"""Duration-weighted fixation maps and image-resampled cross-task similarity."""
from __future__ import annotations

import argparse
import json
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import loadmat
from scipy.ndimage import convolve1d
from skimage.metrics import structural_similarity

from .behavior import load_inclusion, normalize_image_name
from .common import TASKS, canonical_task, fdr_bh


def fixation_probability_map(fixations: pd.DataFrame, width: int = 1920,
                             height: int = 1080, bin_size: int = 4,
                             sigma_bins: float = 6.0) -> np.ndarray:
    """Bin EyeLink events by clipped duration in the [0,4000] ms window.

    X/Y are zero-based screen pixels. Start_ms/End_ms are EyeLink fixation
    event times relative to image onset. The Gaussian kernel has radius
    ceil(3*sigma), and convolution uses zero padding.
    """
    if min(width, height, bin_size) < 1 or sigma_bins <= 0:
        raise ValueError("Positive dimensions, bin size and smoothing sigma required")
    required = ["X", "Y", "Start_ms", "End_ms"]
    if not set(required).issubset(fixations):
        raise ValueError(f"Fixation columns required: {required}")
    values = fixations[required].apply(pd.to_numeric, errors="raise").to_numpy(float)
    if not np.isfinite(values).all():
        raise ValueError("Fixation events must have finite coordinates/times")
    if np.any(values[:, 3] < values[:, 2]):
        raise ValueError("Fixation end precedes start")
    x, y, start, end = values.T
    duration = np.maximum(0, np.minimum(end, 4000) - np.maximum(start, 0))
    keep = (x >= 0) & (x < width) & (y >= 0) & (y < height) & (duration > 0)
    result = np.zeros((int(np.ceil(height / bin_size)), int(np.ceil(width / bin_size))), dtype=float)
    np.add.at(result, (np.floor(y[keep] / bin_size).astype(int),
                      np.floor(x[keep] / bin_size).astype(int)), duration[keep])
    if result.sum() <= 0:
        raise ValueError("No positive-duration in-screen fixation in the 4-s window")
    radius = max(1, int(np.ceil(3 * sigma_bins)))
    support = np.arange(-radius, radius + 1)
    kernel = np.exp(-support**2 / (2 * sigma_bins**2))
    kernel /= kernel.sum()
    result = convolve1d(convolve1d(result, kernel, axis=0, mode="constant", cval=0),
                        kernel, axis=1, mode="constant", cval=0)
    return result / result.sum()


def build_group_maps(trials: pd.DataFrame, fixations: pd.DataFrame, output_dir: Path,
                     inclusion: set | None = None,
                     min_valid_ratio: float = 0.70, min_duration_ms: float = 3900,
                     bin_size: int = 4, sigma_bins: float = 6) -> tuple[pd.DataFrame, pd.DataFrame]:
    """QC each trial, normalize, average within person, then across people.

    Calibration, valid-sample ratio and recording duration determine trial
    eligibility. Participant means receive equal weight in each group map.
    """
    keys = ["Group", "Subject", "Run", "Trial"]
    qc_columns = ["CalibrationOK", "GazeValidRatio", "ValidRecordingDuration_ms"]
    if not set(keys + ["Task", "ImageName"] + qc_columns).issubset(trials):
        raise ValueError(f"Trial CSV requires {keys + ['Task', 'ImageName'] + qc_columns}")
    if not set(keys + ["X", "Y", "Start_ms", "End_ms"]).issubset(fixations):
        raise ValueError("Fixation CSV requires Group,Subject,Run,Trial,X,Y,Start_ms,End_ms")
    trials, fixations = trials.copy(), fixations.copy()
    for table in (trials, fixations):
        table["Group"] = table.Group.astype(str).str.strip()
        table["Subject"] = table.Subject.map(lambda x: f"{int(x):03d}")
        for key in ("Run", "Trial"):
            table[key] = pd.to_numeric(table[key], errors="raise").astype(int)
    if trials.duplicated(keys).any():
        raise ValueError("Duplicate trial IDs")
    event_groups = dict(tuple(fixations.groupby(keys, sort=False)))
    trial_keys = set(map(tuple, trials[keys].to_numpy()))
    if set(event_groups) - trial_keys:
        raise ValueError("Fixation events reference trials absent from the trial metadata")
    # Repeated presentations never give a participant greater group weight.
    subject_sums, subject_counts, qc = {}, {}, []
    group_sums, group_counts = {}, {}

    def flush_subject() -> None:
        # Aggregate one participant at a time to limit peak memory use.
        for subject_key, total in subject_sums.items():
            key = subject_key[:2]
            if key not in group_sums:
                group_sums[key], group_counts[key] = np.zeros_like(total), 0
            group_sums[key] += total / subject_counts[subject_key]
            group_counts[key] += 1
        subject_sums.clear()
        subject_counts.clear()

    last_person = None
    for row in trials.sort_values(["Group", "Subject", "Run", "Trial"]).to_dict("records"):
        trial_key = tuple(row[k] for k in keys)
        if trial_key[:2] != last_person:
            flush_subject()
            last_person = trial_key[:2]
        task = canonical_task(row["Task"])
        image = normalize_image_name(row["ImageName"])
        ratio = float(row["GazeValidRatio"])
        duration = float(row["ValidRecordingDuration_ms"])
        calibration = str(row["CalibrationOK"]).strip().lower()
        if calibration not in ("true", "false", "1", "0", "yes", "no"):
            raise ValueError(f"Invalid CalibrationOK for {trial_key}")
        if not np.isfinite([ratio, duration]).all() or not 0 <= ratio <= 1 or not 0 <= duration <= 4000:
            raise ValueError(f"Missing/out-of-range trial QC for {trial_key}")
        status = "INCLUDED"
        if inclusion is not None and trial_key[:3] not in inclusion and (*trial_key[:2], None) not in inclusion:
            status = "EXCLUDED_BY_MANIFEST"
        elif calibration in ("false", "0", "no"):
            status = "CALIBRATION_FAILED"
        elif ratio < min_valid_ratio:
            status = "LOW_VALID_SAMPLE_RATIO"
        elif duration < min_duration_ms:
            status = "SHORT_VALID_RECORDING"
        elif trial_key not in event_groups:
            status = "NO_FIXATION_EVENTS"
        if status == "INCLUDED":
            try:
                heatmap = fixation_probability_map(event_groups[trial_key], bin_size=bin_size,
                                                  sigma_bins=sigma_bins)
            except ValueError as exc:
                # Trials without a usable fixation receive an explicit exclusion status.
                if "No positive-duration" not in str(exc):
                    raise
                status = "NO_VALID_IN_WINDOW_FIXATION"
            else:
                subject_key = (task, image, *trial_key[:2])
                if subject_key not in subject_sums:
                    subject_sums[subject_key] = np.zeros_like(heatmap)
                    subject_counts[subject_key] = 0
                subject_sums[subject_key] += heatmap
                subject_counts[subject_key] += 1
        qc.append({**{k: row[k] for k in keys}, "Task": task, "ImageName": image,
                   "GazeValidRatio": ratio, "ValidRecordingDuration_ms": duration, "Status": status})
    flush_subject()
    output_dir.mkdir(parents=True, exist_ok=True)
    map_dir = output_dir / "numeric_maps"
    map_dir.mkdir(exist_ok=True)
    index = []
    for (task, image), total in sorted(group_sums.items()):
        relative_path = Path("numeric_maps") / f"{task}_{Path(image).stem}.npy"
        np.save(output_dir / relative_path, total / group_counts[(task, image)])
        index.append(dict(Task=task, ImageName=image, HeatmapFile=relative_path.as_posix(),
                          NSubjects=group_counts[(task, image)]))
    if not index:
        raise ValueError("No trial survived gaze QC; no numeric group maps created")
    return pd.DataFrame(index), pd.DataFrame(qc)


def load_numeric_map(path: Path) -> np.ndarray:
    """Read a nonnegative 2-D numeric map and normalize it to unit mass."""
    if path.suffix.lower() == ".npy":
        heatmap = np.load(path, allow_pickle=False)
    elif path.suffix.lower() == ".npz":
        with np.load(path, allow_pickle=False) as archive:
            heatmap = archive["H"]
    elif path.suffix.lower() == ".mat":
        heatmap = loadmat(path)["H"]
    else:
        raise ValueError(f"Only numeric .npy/.npz/.mat H maps accepted: {path}")
    heatmap = np.asarray(heatmap, dtype=float)
    if heatmap.ndim != 2 or not np.isfinite(heatmap).all() or np.any(heatmap < 0) or heatmap.sum() <= 0:
        raise ValueError(f"Expected finite nonnegative 2-D map with positive mass: {path}")
    return heatmap / heatmap.sum()


def load_group_index(index_path: Path, expected_images: int = 400
                     ) -> dict[tuple[str, str], np.ndarray]:
    """One already aggregated map per Task/ImageName; relative paths use CSV location."""
    index = pd.read_csv(index_path)
    if not {"Task", "ImageName", "HeatmapFile"}.issubset(index):
        raise ValueError("Map index requires Task,ImageName,HeatmapFile")
    means, shape = {}, None
    for row in index.to_dict("records"):
        key = canonical_task(row["Task"]), normalize_image_name(row["ImageName"])
        if key in means:
            raise ValueError(f"Duplicate group map; aggregate trial maps before similarity: {key}")
        path = Path(row["HeatmapFile"])
        path = path if path.is_absolute() else index_path.parent / path
        value = load_numeric_map(path)
        shape = shape or value.shape
        if value.shape != shape:
            raise ValueError(f"Map shape mismatch {value.shape} vs {shape}: {path}")
        means[key] = value
    expected = {f"image{i}.png" for i in range(1, expected_images + 1)}
    for task in TASKS:
        actual = {image for t, image in means if t == task}
        if actual != expected:
            raise ValueError(f"{task}: expected exactly image1..image{expected_images}; "
                             f"missing={sorted(expected - actual)}, extra={sorted(actual - expected)}")
    return means


def r_from_image_statistics(statistics: np.ndarray) -> float:
    """Pearson r from pooled per-image sufficient statistics; no pixel resampling."""
    sx, sy, sxx, syy, sxy, n = statistics.sum(axis=0)
    vx, vy = sxx - sx * sx / n, syy - sy * sy / n
    if min(vx, vy) <= 0:
        return np.nan
    return float(np.clip((sxy - sx * sy / n) / np.sqrt(vx * vy), -1, 1))


def image_ssim(x: np.ndarray, y: np.ndarray) -> float:
    """Compute SSIM with a shared pairwise range and largest odd <=7 window."""
    value_range = max(x.max(), y.max()) - min(x.min(), y.min())
    if value_range <= np.finfo(float).eps:
        return 1.0 if np.allclose(x, y) else np.nan
    win = min(7, *x.shape)
    win -= 1 - win % 2
    if win < 3:
        raise ValueError("SSIM needs map dimensions >=3")
    return float(structural_similarity(x, y, data_range=value_range, win_size=win))


def task_similarity(means: dict, n_boot: int = 5000,
                    seed: int = 20260717) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Six comparisons; CI and one-sided empirical P use whole-image bootstrap.

    P=(# bootstrap estimates <=0 +1)/(B+1) quantifies the empirical tail
    probability for positive similarity. BH correction is separate for
    Pearson correlation and SSIM.
    """
    if n_boot < 1:
        raise ValueError("n_boot must be positive")
    image_sets = [{image for task, image in means if task == t} for t in TASKS]
    if not image_sets[0] or any(s != image_sets[0] for s in image_sets):
        raise ValueError("All four tasks must have the same common image set")
    images = sorted(image_sets[0], key=lambda name: int(Path(name).stem[5:]))
    if len(images) < 3:
        raise ValueError("At least three common images required")
    results, per_image = [], []
    # Use deterministic image-resampling streams for each similarity metric.
    rng_r, rng_s = np.random.default_rng(seed), np.random.default_rng(seed)
    for task1, task2 in combinations(TASKS, 2):
        statistics, ssim_values = [], []
        for image in images:
            x, y = means[(task1, image)], means[(task2, image)]
            if x.shape != y.shape:
                raise ValueError("Paired maps differ in shape")
            statistics.append([x.sum(), y.sum(), np.square(x).sum(), np.square(y).sum(),
                               (x * y).sum(), x.size])
            score = image_ssim(x, y)
            ssim_values.append(score)
            per_image.append(dict(Task1=task1, Task2=task2, ImageName=image, SSIM=score))
        statistics, ssim_values = np.asarray(statistics), np.asarray(ssim_values)
        if not np.isfinite(ssim_values).all():
            raise ValueError(f"Nonfinite SSIM for {task1}/{task2}; no image silently discarded")
        boot_r, boot_s = np.empty(n_boot), np.empty(n_boot)
        for b in range(n_boot):
            boot_r[b] = r_from_image_statistics(statistics[rng_r.integers(len(images), size=len(images))])
            boot_s[b] = ssim_values[rng_s.integers(len(images), size=len(images))].mean()
        for metric, estimate, boot in (("Pearson", r_from_image_statistics(statistics), boot_r),
                                       ("SSIM", ssim_values.mean(), boot_s)):
            valid_boot = boot[np.isfinite(boot)]
            if not np.isfinite(estimate) or len(valid_boot) != n_boot:
                raise ValueError(f"Degenerate {metric} bootstrap for {task1}/{task2}")
            low, high = np.percentile(valid_boot, [2.5, 97.5])
            results.append(dict(Task1=task1, Task2=task2, Metric=metric, Estimate=float(estimate),
                                CI95_low=float(low), CI95_high=float(high),
                                P_bootstrap_greater_than_zero=float((np.sum(valid_boot <= 0) + 1) / (n_boot + 1)),
                                N_Images=len(images), N_Bootstrap=n_boot, Seed=seed))
    result = pd.DataFrame(results)
    result["Q_FDR_6Pairs_within_metric"] = np.nan
    for metric in ("Pearson", "SSIM"):
        selected = result.Metric == metric
        result.loc[selected, "Q_FDR_6Pairs_within_metric"] = fdr_bh(
            result.loc[selected, "P_bootstrap_greater_than_zero"].to_numpy()
        )
    return result, pd.DataFrame(per_image)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    build = commands.add_parser("from-fixations", help="Build numeric maps from aligned EyeLink fixation events")
    build.add_argument("--trials", type=Path, required=True)
    build.add_argument("--fixations", type=Path, required=True)
    build.add_argument("--participants", type=Path, required=True)
    build.add_argument("--min-valid-ratio", type=float, default=0.70)
    build.add_argument("--min-duration-ms", type=float, default=3900)
    build.add_argument("--bin-size", type=int, default=4)
    build.add_argument("--sigma-bins", type=float, default=6)
    similarity = commands.add_parser("similarity", help="Analyze already aggregated, numeric probability maps")
    similarity.add_argument("--index", type=Path, required=True)
    similarity.add_argument("--n-bootstrap", type=int, default=5000)
    similarity.add_argument("--seed", type=int, default=20260717)
    similarity.add_argument("--expected-images", type=int, default=400)
    for command in (build, similarity):
        command.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    metadata = {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()}
    if args.command == "from-fixations":
        if not 0 <= args.min_valid_ratio <= 1 or not 0 <= args.min_duration_ms <= 4000:
            raise ValueError("Gaze QC thresholds outside valid range")
        index, qc = build_group_maps(pd.read_csv(args.trials, dtype={"Subject": str}),
                                    pd.read_csv(args.fixations, dtype={"Subject": str}), args.output_dir,
                                    load_inclusion(args.participants),
                                    args.min_valid_ratio, args.min_duration_ms, args.bin_size, args.sigma_bins)
        index.to_csv(args.output_dir / "group_map_index.csv", index=False)
        qc.to_csv(args.output_dir / "gaze_trial_qc.csv", index=False)
        metadata["status_counts"] = {k: int(v) for k, v in qc.Status.value_counts().items()}
    else:
        means = load_group_index(args.index, args.expected_images)
        result, per_image = task_similarity(means, args.n_bootstrap, args.seed)
        result.to_csv(args.output_dir / "gaze_task_similarity.csv", index=False)
        per_image.to_csv(args.output_dir / "image_level_ssim.csv", index=False)
        metadata["map_level"] = "task_by_image_group_mean"
    (args.output_dir / "analysis_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
