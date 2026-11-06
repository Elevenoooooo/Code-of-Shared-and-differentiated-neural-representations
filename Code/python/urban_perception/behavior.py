"""Image-level behavioural profiles and cross-task rating correlations."""
from __future__ import annotations

import argparse
import json
import re
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

from .common import TASKS, canonical_task, fdr_bh

GROUP_RUN_TASK = {
    "Group1": ("Safety", "Vitality", "Beauty", "Imageability"),
    "Group2": ("Vitality", "Safety", "Imageability", "Beauty"),
    "Group3": ("Beauty", "Imageability", "Safety", "Vitality"),
    "Group4": ("Imageability", "Beauty", "Vitality", "Safety"),
}
FILENAME = re.compile(
    r"^(Group[1-4])_(\d+)_run([1-4])_([a-z]+)_trials\.csv$", re.I
)


def basename(value: object) -> str:
    """Read Windows or POSIX paths independent of the host operating system."""
    return str(value).strip().replace("\\", "/").rsplit("/", 1)[-1].lower()


def normalize_image_name(value: object) -> str:
    """Extract and validate the imageN.png identifier from a stimulus path."""
    name = basename(value)
    if not re.fullmatch(r"image[1-9]\d*\.png", name):
        raise ValueError(f"Expected an imageN.png identifier: {value!r}")
    return name


def load_inclusion(path: str | Path) -> set[tuple[str, str, int | None]]:
    """CSV: Group,Subject and optional Run,Include (true/false or 1/0).

    A missing Run applies the decision to all runs. Subject numbers are local to
    Group, so Group1_001 and Group2_001 are separate participants.
    """
    frame = pd.read_csv(path, dtype=str, encoding="utf-8-sig").fillna("")
    if not {"Group", "Subject"}.issubset(frame):
        raise ValueError("Participant CSV requires Group and Subject")
    result = set()
    seen = set()
    for row in frame.to_dict("records"):
        group = str(row["Group"]).strip()
        if group not in GROUP_RUN_TASK:
            raise ValueError(f"Unknown Group: {group}")
        subject_text = row["Subject"].strip()
        if not subject_text.isdigit():
            raise ValueError(f"Subject must be a within-group number: {subject_text}")
        subject = f"{int(subject_text):03d}"
        run_text = row.get("Run", "").strip().lower().removeprefix("run")
        run = int(run_text) if run_text else None
        if run is not None and run not in (1, 2, 3, 4):
            raise ValueError(f"Invalid Run: {run}")
        key = (group, subject, run)
        if key in seen:
            raise ValueError(f"Duplicate participant/run decision: {key}")
        seen.add(key)
        include = row.get("Include", "true").strip().lower()
        if include not in ("true", "false", "1", "0", "yes", "no"):
            raise ValueError(f"Invalid Include: {include}")
        if include in ("true", "1", "yes"):
            result.add(key)
    if not result:
        raise ValueError("No included participant/run entries")
    # Mixed participant-wide and run-specific decisions would be ambiguous.
    for group, subject, run in seen:
        if run is not None and (group, subject, None) in seen:
            raise ValueError(f"Use either participant-wide OR run decisions: {group}_{subject}")
    return result


def collect_binary_ratings(behavior_dir: str | Path, inclusion: set | None = None,
                           expected_images: int = 400) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read formal trials and recode key 1=-1/key 2=+1.

    Path_raw identifies the presented image. Numeric thisTrialN_raw selects
    formal trials, excluding PsychoPy footers. Invalid keys are retained as
    NaN in the trial table and counted in the run-level quality summary.
    """
    root = Path(behavior_dir)
    paths = sorted(root.rglob("*_trials.csv"))
    if not paths:
        raise FileNotFoundError(f"No *_trials.csv under {root}")
    trials, quality, found = [], [], set()
    canonical = {f"image{i}.png" for i in range(1, expected_images + 1)}
    for path in paths:
        match = FILENAME.fullmatch(path.name)
        if not match:
            raise ValueError(f"Unrecognized trial filename: {path.name}")
        group = f"Group{match[1][-1]}"
        subject, run = f"{int(match[2]):03d}", int(match[3])
        task = canonical_task(match[4])
        if GROUP_RUN_TASK[group][run - 1] != task:
            raise ValueError(f"Filename task disagrees with Latin-square schedule: {path}")
        key = (group, subject, run)
        if key in found:
            raise ValueError(f"Duplicate participant/run file: {key}")
        found.add(key)
        included = inclusion is None or key in inclusion or (group, subject, None) in inclusion
        if not included:
            quality.append(dict(Group=group, Subject=subject, Run=run, Task=task,
                                File=str(path.relative_to(root)), Status="EXCLUDED_BY_MANIFEST"))
            continue
        raw = pd.read_csv(path, dtype=str, encoding="utf-8-sig")
        raw.columns = raw.columns.str.strip()
        required = {"Path_raw", "key_resp_sc.keys_raw", "thisTrialN_raw"}
        if not required.issubset(raw):
            raise ValueError(f"Missing columns {required - set(raw)} in {path}")
        formal_number = pd.to_numeric(raw.thisTrialN_raw, errors="coerce")
        formal = raw.loc[formal_number.notna()].copy()
        formal["Trial"] = formal_number[formal_number.notna()].astype(int) + 1
        if formal["Trial"].duplicated().any() or len(formal) > 100:
            raise ValueError(f"Duplicate or >100 formal trials in {path}")
        formal["ImageName"] = formal.Path_raw.map(normalize_image_name)
        if not formal.ImageName.isin(canonical).all():
            raise ValueError(f"Image outside image1..image{expected_images}: {path}")
        formal["RawKey"] = pd.to_numeric(formal["key_resp_sc.keys_raw"], errors="coerce")
        formal["SignedRating"] = formal.RawKey.map({1: -1.0, 2: 1.0})
        formal["ResponseTime_s"] = pd.to_numeric(
            formal.get("key_resp_sc.rt_raw", pd.Series(np.nan, index=formal.index)), errors="coerce"
        )
        out = formal[["Trial", "ImageName", "RawKey", "SignedRating", "ResponseTime_s"]].copy()
        out.insert(0, "Task", task)
        out.insert(0, "Run", run)
        out.insert(0, "Subject", subject)
        out.insert(0, "Group", group)
        trials.append(out)
        quality.append(dict(Group=group, Subject=subject, Run=run, Task=task,
                            File=str(path.relative_to(root)), Status="INCLUDED",
                            NFormalTrials=len(formal), NValidKeys=int(formal.SignedRating.notna().sum()),
                            NInvalidKeys=int(formal.SignedRating.isna().sum()),
                            NRepeatedImages=int(formal.ImageName.duplicated().sum())))
    if inclusion is not None:
        for group, subject, run in sorted(inclusion, key=str):
            for selected_run in ((1, 2, 3, 4) if run is None else (run,)):
                if (group, subject, selected_run) not in found:
                    raise FileNotFoundError(f"Included run has no data: {group}_{subject}, run{selected_run}")
    if not trials:
        raise ValueError("No included behavioural files")
    return pd.concat(trials, ignore_index=True), pd.DataFrame(quality)


def build_image_profiles(trials: pd.DataFrame) -> pd.DataFrame:
    """Equal participant weighting: trial repeats first average within person."""
    valid = trials.dropna(subset=["SignedRating"])
    subject_means = valid.groupby(["Group", "Subject", "Task", "ImageName"], as_index=False).agg(
        SignedRating=("SignedRating", "mean")
    )
    profile = subject_means.groupby(["Task", "ImageName"], as_index=False).agg(
        NetRatingMean=("SignedRating", "mean"), NetRatingSD=("SignedRating", "std"),
        NSubjects=("SignedRating", "count")
    )
    profile["PositiveRate"] = (profile.NetRatingMean + 1.0) / 2.0
    profile["NetScore_100"] = 100 * profile.NetRatingMean
    return profile


def rating_correlations(wide: pd.DataFrame) -> pd.DataFrame:
    """Pearson correlations across images; BH across the six task pairs."""
    if not set(TASKS).issubset(wide.columns):
        raise ValueError(f"Expected task columns: {TASKS}")
    rows = []
    for task1, task2 in combinations(TASKS, 2):
        pair = wide[[task1, task2]].replace([np.inf, -np.inf], np.nan).dropna()
        if len(pair) >= 3 and pair[task1].std() > 0 and pair[task2].std() > 0:
            r, p = pearsonr(pair[task1], pair[task2])
        else:
            r, p = np.nan, np.nan
        rows.append(dict(Task1=task1, Task2=task2, Pearson_r=r,
                         P_uncorrected=p, N_Images=len(pair)))
    result = pd.DataFrame(rows)
    result["Q_FDR_6Pairs"] = fdr_bh(result.P_uncorrected.to_numpy())
    return result


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    binary = commands.add_parser("binary", help="Summarize fMRI behavioural responses")
    binary.add_argument("--behavior-dir", type=Path, required=True)
    binary.add_argument("--participants", type=Path, required=True, help="Participant inclusion CSV")
    binary.add_argument("--expected-images", type=int, default=400)
    online = commands.add_parser("online", help="Correlate an existing image-level online rating table")
    online.add_argument("--ratings", type=Path, required=True)
    online.add_argument("--image-column", default="File_Name")
    for command in (binary, online):
        command.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.command == "binary":
        inclusion = load_inclusion(args.participants)
        trials, quality = collect_binary_ratings(args.behavior_dir, inclusion, args.expected_images)
        profile = build_image_profiles(trials)
        wide = profile.pivot(index="ImageName", columns="Task", values="PositiveRate").reindex(columns=TASKS)
        trials.to_csv(args.output_dir / "binary_rating_trials.csv", index=False)
        quality.to_csv(args.output_dir / "behavior_qc.csv", index=False)
        profile.to_csv(args.output_dir / "image_level_rating_profiles.csv", index=False)
        metadata = dict(analysis="binary_ratings",
                        n_unique_participants=int(trials[["Group", "Subject"]].drop_duplicates().shape[0]),
                        n_subject_runs=int(trials[["Group", "Subject", "Run"]].drop_duplicates().shape[0]),
                        n_images=int(trials.ImageName.nunique()), image_column="Path_raw",
                        included_manifest=str(args.participants))
    else:
        raw = pd.read_csv(args.ratings)
        raw = raw.rename(columns={"safety": "Safety", "liveliness": "Vitality",
                                  "beauty": "Beauty", "imageability": "Imageability"})
        if args.image_column not in raw:
            raise ValueError("Online table requires a unique image ID column")
        raw[args.image_column] = raw[args.image_column].map(normalize_image_name)
        if raw[args.image_column].duplicated().any():
            raise ValueError("Online table requires a unique image ID column")
        wide = raw.set_index(args.image_column).reindex(columns=TASKS).apply(pd.to_numeric, errors="raise")
        if wide.isna().all().any():
            raise ValueError("Online table requires ratings for all four tasks")
        metadata = dict(analysis="online_image_means", input=str(args.ratings), n_images=len(wide),
                        aggregation_level="image_mean")
    rating_correlations(wide).to_csv(args.output_dir / "rating_correlations.csv", index=False)
    wide.to_csv(args.output_dir / "rating_profiles_wide.csv")
    (args.output_dir / "analysis_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
