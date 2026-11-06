"""Map-based activation summaries, profile similarity, and SWS reconstruction.

Uses unthresholded statistics and separately supplied q<.05 masks. NIfTI
affines must match; no resampling is performed silently. These functions
implement whole-ROI signed-T summaries and zero-sum task coordinates.
"""

import argparse
from itertools import combinations
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from .common import TASKS, TASK_PAIRS, fdr_bh


def sws_from_voxels(signed_t, significant):
    """Coverage times adjusted Fisher-Pearson skewness of ALL finite ROI T's.

    The supplied mask must come from one two-sided A-B test, not two separately
    corrected one-sided maps. This descriptive index is not a direction test.
    """
    t = np.asarray(signed_t, dtype=float)
    mask = np.asarray(significant)
    if t.shape != mask.shape or not np.all(np.isfinite(mask)):
        raise ValueError("T and significance arrays must align; mask must be finite.")
    values = t[np.isfinite(t)]
    if values.size < 3 or np.ptp(values) == 0:
        raise ValueError("SWS requires at least three finite nonconstant T values.")
    coverage = np.count_nonzero(mask[np.isfinite(t)] > 0) / values.size
    skewness = float(stats.skew(values, bias=False))
    return {"N_Voxels": values.size, "Coverage": coverage,
            "AdjustedSkewness": skewness, "SWS": coverage * skewness}


def reconstruct_coordinates(pairwise, rois=("VVC", "PPC", "vmPFC")):
    """Fit six oriented pairwise SWS differences per ROI under sum(x)=0.

    Required columns: TaskA, TaskB, ROI, SWS. Every unordered task pair occurs
    exactly once per ROI. Reversed input contrasts are sign-corrected.
    Returns coordinates and pairwise residuals; no PCA or axis standardization.
    """
    required = {"TaskA", "TaskB", "ROI", "SWS"}
    if not required.issubset(pairwise.columns):
        raise ValueError(f"Missing columns: {required - set(pairwise.columns)}")
    design = np.zeros((6, 4))
    for row, (a, b) in enumerate(TASK_PAIRS):
        design[row, TASKS.index(a)] = 1
        design[row, TASKS.index(b)] = -1
    coordinates = pd.DataFrame({"Task": TASKS})
    residual_rows = []
    for roi in rois:
        subset = pairwise.loc[pairwise.ROI.eq(roi)]
        observed = {}
        for row in subset.itertuples(index=False):
            a, b = row.TaskA, row.TaskB
            if a not in TASKS or b not in TASKS or a == b or not np.isfinite(row.SWS):
                raise ValueError(f"Invalid pair or SWS for {roi}")
            key = (a, b) if TASKS.index(a) < TASKS.index(b) else (b, a)
            if key in observed:
                raise ValueError(f"Duplicate task pair {key} for {roi}")
            observed[key] = float(row.SWS) if key == (a, b) else -float(row.SWS)
        if set(observed) != set(TASK_PAIRS):
            raise ValueError(f"Need all six task pairs for {roi}")
        y = np.array([observed[pair] for pair in TASK_PAIRS])
        # The graph incidence matrix has rank 3; minimum norm fixes zero sum.
        x = np.linalg.lstsq(design, y, rcond=None)[0]
        x -= x.mean()
        coordinates[roi] = x
        fitted = design @ x
        for (a, b), obs, fit in zip(TASK_PAIRS, y, fitted):
            residual_rows.append({"ROI": roi, "TaskA": a, "TaskB": b,
                                  "Observed": obs, "Fitted": fit, "Residual": obs - fit})
    return coordinates, pd.DataFrame(residual_rows)


def coordinate_distances(coordinates, rois=("VVC", "PPC", "vmPFC")):
    """Euclidean distances in the three ROI SWS units."""
    values = coordinates.set_index("Task").loc[list(TASKS), list(rois)]
    return pd.DataFrame([
        {"TaskA": a, "TaskB": b, "Distance": float(np.linalg.norm(values.loc[a] - values.loc[b]))}
        for a, b in TASK_PAIRS
    ])


def profile_similarity(profiles):
    """Pearson similarity of four aligned profiles, BH across six pairs.

    Rows index identically ordered voxels or 105 unique FC edges. Missing or
    constant profiles are errors, rather than changing the sample per pair.
    Pearson P values assume independent observations; spatial/network
    dependence applies to the parametric interpretation.
    """
    array = profiles.loc[:, list(TASKS)].to_numpy(dtype=float)
    if array.shape[0] < 3 or not np.isfinite(array).all() or np.any(np.ptp(array, axis=0) == 0):
        raise ValueError("Profiles must contain >=3 finite, nonconstant aligned observations.")
    rows = []
    for a, b in TASK_PAIRS:
        r, p = stats.pearsonr(array[:, TASKS.index(a)], array[:, TASKS.index(b)])
        rows.append({"TaskA": a, "TaskB": b, "r": r, "p": p, "N": array.shape[0]})
    table = pd.DataFrame(rows)
    table["q"] = fdr_bh(table.p)
    return table


def load_aligned_images(paths):
    """Read same-grid NIfTIs; nibabel is an optional neuroimaging dependency."""
    import nibabel as nib
    images = [nib.load(str(p)) for p in paths]
    if not images:
        raise ValueError("No NIfTI files supplied")
    reference = images[0]
    if len(reference.shape) != 3:
        raise ValueError("Expected 3D statistical images")
    for image in images[1:]:
        if image.shape != reference.shape or not np.allclose(image.affine, reference.affine, rtol=0, atol=1e-5):
            raise ValueError("NIfTI grids differ; explicitly reslice atlas using nearest neighbour.")
    return [im.get_fdata(dtype=np.float64) for im in images]


def sws_from_maps(map_dir, atlas_path, roi_table):
    """Extract SWS from signed A-B maps and their two-sided FDR masks."""
    rows = []
    for a, b in TASK_PAIRS:
        stem = f"Stim_gt_Fixation_{a}_gt_{b}"
        atlas, t, mask = load_aligned_images([
            atlas_path, Path(map_dir) / f"{stem}_statmap.nii",
            Path(map_dir) / f"{stem}_sigmask_q050.nii"])
        for roi in roi_table.itertuples(index=False):
            selected = atlas == roi.ROI_ID
            result = sws_from_voxels(t[selected], mask[selected])
            rows.append({"TaskA": a, "TaskB": b, "ROI": roi.ROI, **result})
    return pd.DataFrame(rows)


def activation_from_maps(map_dir, atlas_path, roi_table):
    """Summarize retained ROIs using supplied unthresholded maps and q masks.

    Outputs: ROI differentiation, common activation counts, concatenated-profile
    similarity (six tests), and ROI-specific profile similarity (90 tests).
    Input maps represent Stimulus > Fixation and the task main effect.
    """
    folder = Path(map_dir)
    paths = [atlas_path]
    paths += [folder / f"Stim_gt_Fixation_OneSample_{t}_gt_0_statmap.nii" for t in TASKS]
    paths += [folder / f"Stim_gt_Fixation_OneSample_{t}_gt_0_sigmask_q050.nii" for t in TASKS]
    paths += [folder / "Stim_gt_Fixation_TaskMainEffect_F_statmap.nii",
              folder / "Stim_gt_Fixation_TaskMainEffect_F_sigmask_q050.nii"]
    atlas, *arrays = load_aligned_images(paths)
    t_maps, masks = arrays[:4], arrays[4:8]
    f_map, f_sig = arrays[8:]
    if any(not np.isfinite(m).all() for m in [*masks, f_sig]):
        raise ValueError("Significance masks contain missing values")
    common = np.logical_and.reduce([mask > 0 for mask in masks])
    summaries, profiles, roi_correlations = [], [], []
    for roi in roi_table.itertuples(index=False):
        selected = atlas == roi.ROI_ID
        if not selected.any() or not np.isfinite(f_map[selected]).all():
            raise ValueError(f"Missing or invalid statistical values for {roi.ROI}")
        frame = pd.DataFrame({task: values[selected] for task, values in zip(TASKS, t_maps)})
        comparisons = profile_similarity(frame)
        comparisons.insert(0, "ROI", roi.ROI)
        roi_correlations.append(comparisons)
        profiles.append(frame)
        summaries.append({"ROI": roi.ROI, "N_Voxels": int(selected.sum()),
                          "Mean_F": float(f_map[selected].mean()),
                          "TaskEffect_FDR_Coverage": float((f_sig[selected] > 0).mean()),
                          "Conjunction_N_Voxels": int(common[selected].sum()),
                          "Conjunction_Coverage": float(common[selected].mean())})
    by_roi = pd.concat(roi_correlations, ignore_index=True)
    by_roi["q_across_all_ROI_pairs"] = fdr_bh(by_roi.p)
    return pd.DataFrame(summaries), profile_similarity(pd.concat(profiles)), by_roi


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    sws = commands.add_parser("sws", help="Reconstruct Fig. 4 SWS coordinates from a long CSV")
    sws.add_argument("--pairs", required=True, type=Path)
    sws.add_argument("--out", required=True, type=Path)
    similarity = commands.add_parser("similarity", help="Correlate aligned voxel or FC edge profiles")
    similarity.add_argument("--profiles", required=True, type=Path)
    similarity.add_argument("--out", required=True, type=Path)
    maps = commands.add_parser("activation-maps", help="Summarize supplied group maps with a same-grid merged atlas")
    maps.add_argument("--maps", required=True, type=Path)
    maps.add_argument("--atlas", required=True, type=Path)
    maps.add_argument("--rois", required=True, type=Path)
    maps.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    if args.command == "sws":
        coordinates, residuals = reconstruct_coordinates(pd.read_csv(args.pairs))
        coordinates.to_csv(args.out / "sws_coordinates.csv", index=False)
        residuals.to_csv(args.out / "sws_reconstruction_residuals.csv", index=False)
        coordinate_distances(coordinates).to_csv(args.out / "sws_distances.csv", index=False)
    elif args.command == "similarity":
        profile_similarity(pd.read_csv(args.profiles)).to_csv(args.out / "profile_similarity.csv", index=False)
    else:
        rois = pd.read_csv(args.rois)
        rois = rois.loc[rois.Retained15.astype(str).str.casefold().isin(["true", "1"])]
        summary, union, by_roi = activation_from_maps(args.maps, args.atlas, rois)
        summary.to_csv(args.out / "roi_activation_summary.csv", index=False)
        union.to_csv(args.out / "activation_profile_similarity.csv", index=False)
        by_roi.to_csv(args.out / "roi_profile_similarity.csv", index=False)


if __name__ == "__main__":
    main()
