"""Figure 5: environmental-feature RDMs and run-level RSA algorithms.

The run-level neural and rating RDMs contain 99 retained stimulus trials.
Environmental-feature RDMs represent the full set of 400 stimuli and are
subset to each run by their image labels.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import re
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


# Ten simultaneous environmental predictors in the Figure 5 RSA-GLM.
FEATURE_FILES = {
    "Colorfulness": "RDM_colorfulness.csv",
    "Visual Complexity": "RDM_visual_complexity.csv",
    "Spatial Layering": "RDM_spatial_layering.csv",
    "Green View Index": "RDM_green_view_index.csv",
    "Sky Visibility": "RDM_sky_visibility.csv",
    "Enclosure": "RDM_enclosure.csv",
    "Crowd Attraction": "RDM_crowd_attraction.csv",
    "Vehicle Interference": "RDM_vehicle_interference.csv",
    "Street Facility Convenience": "RDM_street_facility_convenience.csv",
    "Urbanization Intensity Index": "RDM_urbanization_intensity_index.csv",
}


def _stem(value: object) -> str:
    return Path(str(value).strip().replace("\\", "/")).stem.lower()


def image_ids(values) -> list[str]:
    """Normalize standard image1.png through image400.png stimulus labels."""
    keys = [_stem(x) for x in values]
    if any(re.fullmatch(r"image(?:[1-9]\d*)", x) is None for x in keys):
        raise ValueError("Stimulus labels must use image1.png through image400.png")
    if any(not 1 <= int(x[5:]) <= 400 for x in keys):
        raise ValueError("Image ID outside the 400-image stimulus set")
    return keys


def image_order(image_manifest: str | Path | None = None) -> list[str]:
    """Return complete stimulus labels, optionally in experiment-index order."""
    if image_manifest is None:
        return [f"image{i}" for i in range(1, 401)]
    table = pd.read_csv(image_manifest, encoding="utf-8-sig")
    if not {"File_Name", "Img_ID"}.issubset(table.columns):
        raise ValueError("Image manifest requires File_Name and Img_ID")
    indices = pd.to_numeric(table.Img_ID, errors="raise").to_numpy(float)
    if not np.isfinite(indices).all() or len(set(indices)) != len(indices):
        raise ValueError("Img_ID must contain finite unique experiment indices")
    order = image_ids(table.iloc[np.argsort(indices)].File_Name)
    if len(order) != 400 or set(order) != {f"image{i}" for i in range(1, 401)}:
        raise ValueError("Image manifest must contain every stimulus exactly once")
    return order


def validate_rdm(rdm: np.ndarray, n: int | None = None) -> np.ndarray:
    rdm = np.asarray(rdm, dtype=float)
    if rdm.ndim != 2 or rdm.shape[0] != rdm.shape[1]:
        raise ValueError("RDM must be square")
    if n is not None and rdm.shape != (n, n):
        raise ValueError(f"Expected {n} x {n} RDM, found {rdm.shape}")
    if not np.isfinite(rdm).all():
        raise ValueError("Objective-feature RDM contains nonfinite values")
    if not np.allclose(rdm, rdm.T, atol=1e-8, rtol=0):
        raise ValueError("RDM is not symmetric")
    if not np.allclose(np.diag(rdm), 0, atol=1e-8, rtol=0):
        raise ValueError("RDM diagonal is not zero")
    if np.any(rdm < -1e-8):
        raise ValueError("Objective-feature distance cannot be negative")
    return rdm


def scalar_rdm(values, *, standardize: bool = False) -> np.ndarray:
    """Absolute scalar differences, with optional sample-SD standardization."""
    values = np.asarray(values, dtype=float).reshape(-1)
    if not np.isfinite(values).all():
        raise ValueError("Missing feature/rating values must be resolved explicitly")
    if standardize:
        sd = values.std(ddof=1)
        if not np.isfinite(sd) or sd <= np.finfo(float).eps:
            raise ValueError("Cannot standardize a constant feature")
        values = (values - values.mean()) / sd
    return np.abs(values[:, None] - values[None, :])


def rdm_vector(rdm: np.ndarray) -> np.ndarray:
    """Strict lower triangle in MATLAB column-major order (no diagonal)."""
    matrix = np.asarray(rdm)
    mask = np.tril(np.ones(matrix.shape, dtype=bool), -1)
    return matrix.ravel(order="F")[mask.ravel(order="F")]


def _zscore_columns(matrix: np.ndarray) -> np.ndarray:
    matrix = np.asarray(matrix, dtype=float)
    sd = matrix.std(axis=0, ddof=1)
    sd = np.where(sd > 0, sd, 1.0)
    return (matrix - matrix.mean(axis=0)) / sd


def compute_vif(predictors: np.ndarray) -> np.ndarray:
    """OLS variance inflation factors, including an intercept in each fit."""
    x = np.asarray(predictors, dtype=float)
    result = np.ones(x.shape[1])
    for j in range(x.shape[1]):
        others = np.delete(x, j, axis=1)
        if others.shape[1] == 0:
            continue
        design = np.column_stack([np.ones(len(x)), others])
        residual = x[:, j] - design @ np.linalg.lstsq(design, x[:, j], rcond=None)[0]
        sst = np.sum((x[:, j] - x[:, j].mean()) ** 2)
        fraction = np.sum(residual ** 2) / sst if sst > 0 else 0
        # A finite epsilon cap handles a perfectly collinear design.
        result[j] = 1 / max(fraction, np.finfo(float).eps)
    return result


def neural_rdm(patterns: np.ndarray, *, min_trials: int = 20,
               min_voxels: int = 10, expected_trials: int = 99,
               require_expected_trials: bool = True) -> tuple[np.ndarray, np.ndarray, int]:
    """Compute 1 - Pearson correlation between trial-wise voxel patterns.

    Rows are trials and columns are voxels. Finite, nonconstant voxels enter
    the correlation calculation. The returned trial mask defines the same
    image subset for every predictor RDM.
    """
    raw = np.asarray(patterns, dtype=float)
    if raw.ndim != 2:
        raise ValueError("Patterns must be trial x voxel")
    if require_expected_trials and raw.shape[0] != expected_trials:
        raise ValueError(f"Expected {expected_trials} retained trials")
    keep = np.isfinite(raw).any(axis=1)
    if keep.sum() < min_trials:
        raise ValueError("Too few valid trials")
    if require_expected_trials and keep.sum() != expected_trials:
        raise ValueError(f"Expected {expected_trials} finite trial patterns")
    x = raw[keep]
    good = np.isfinite(x).all(axis=0) & (np.nanstd(x, axis=0, ddof=1) > np.finfo(float).eps)
    if good.sum() < min_voxels:
        raise ValueError("Too few nonconstant finite voxels")
    x = x[:, good]
    with np.errstate(invalid="ignore", divide="ignore"):
        distance = 1 - np.corrcoef(x)
    distance = (distance + distance.T) / 2
    np.fill_diagonal(distance, 0)
    return distance, keep, int(good.sum())


@dataclass
class GLMResult:
    beta: np.ndarray
    vif: np.ndarray
    unique_r2: np.ndarray
    predictor_kept: np.ndarray
    r2: float
    adjusted_r2: float
    n_pairs: int


def fit_rsa_glm(neural: np.ndarray, predictors: list[np.ndarray], *,
                min_pairs: int = 1000) -> GLMResult:
    """Simultaneous standardized OLS; no pairwise p values enter inference."""
    y = rdm_vector(neural)
    x = np.column_stack([rdm_vector(p) for p in predictors])
    if any(np.shape(p) != np.shape(neural) for p in predictors):
        raise ValueError("All RDMs must use identical trial order and shape")
    finite = np.isfinite(y) & np.isfinite(x).all(axis=1)
    x, y = x[finite], y[finite]
    if len(y) < min_pairs:
        raise ValueError("Too few valid RDM pairs")
    keep = np.std(x, axis=0, ddof=1) > np.finfo(float).eps
    if not keep.any() or y.std(ddof=1) <= np.finfo(float).eps:
        raise ValueError("Constant response or all constant predictors")
    x = _zscore_columns(x[:, keep])
    y = _zscore_columns(y[:, None]).ravel()
    design = np.column_stack([np.ones(len(y)), x])
    beta = np.linalg.lstsq(design, y, rcond=None)[0]
    residual = y - design @ beta
    sst = np.sum((y - y.mean()) ** 2)
    r2 = float(1 - np.sum(residual ** 2) / sst)
    unique = np.empty(x.shape[1])
    for j in range(x.shape[1]):
        reduced = np.delete(design, j + 1, axis=1)
        fitted = reduced @ np.linalg.lstsq(reduced, y, rcond=None)[0]
        unique[j] = r2 - (1 - np.sum((y - fitted) ** 2) / sst)
    adj = 1 - (1 - r2) * (len(y) - 1) / max(len(y) - x.shape[1] - 1, 1)
    return GLMResult(beta[1:], compute_vif(x), unique, keep, r2, float(adj), len(y))


def rating_rsa(neural: np.ndarray, ratings, *, min_pairs: int = 1000) -> tuple[float, float, int]:
    """Spearman rho and clipped Fisher z; inference is across subjects."""
    values = np.asarray(ratings, dtype=float).reshape(-1)
    if np.shape(neural) != (len(values), len(values)):
        raise ValueError("Ratings must follow the exact neural-RDM trial order")
    # Exclude pairs involving missing responses.
    rating_distances = np.abs(values[:, None] - values[None, :])
    y, x = rdm_vector(neural), rdm_vector(rating_distances)
    valid = np.isfinite(y) & np.isfinite(x)
    if valid.sum() < min_pairs:
        raise ValueError("Too few valid RDM pairs")
    rho = float(spearmanr(x[valid], y[valid]).statistic)
    return rho, float(np.arctanh(np.clip(rho, -0.999999, 0.999999))), int(valid.sum())


def read_labeled_rdm(path: str | Path, order: list[str]) -> np.ndarray:
    """Align both RDM axes to the requested standard image-label order."""
    table = pd.read_csv(path, encoding="utf-8-sig")
    rows = image_ids(table.iloc[:, 0])
    cols = image_ids(table.columns[1:])
    if len(set(rows)) != len(rows) or len(set(cols)) != len(cols):
        raise ValueError(f"Duplicate RDM labels: {path}")
    if set(rows) != set(order) or set(cols) != set(order):
        raise ValueError(f"RDM labels do not match the full image set: {path}")
    values = pd.DataFrame(table.iloc[:, 1:].to_numpy(float), index=rows, columns=cols)
    return validate_rdm(values.loc[order, order].to_numpy(), len(order))


def prepare_feature_rdms(*, output_dir: str | Path,
                         feature_csv: str | Path | None = None,
                         rdm_dir: str | Path | None = None,
                         image_manifest: str | Path | None = None) -> None:
    """Construct or align ten 400-image environmental-feature RDMs."""
    if (feature_csv is None) == (rdm_dir is None):
        raise ValueError("Specify exactly one of feature_csv and rdm_dir")
    order = image_order(image_manifest)
    rdms = {}
    aligned = None
    if feature_csv is not None:
        table = pd.read_csv(feature_csv, encoding="utf-8-sig")
        idcol = next((c for c in ("ImageName", "File_Name") if c in table), None)
        if idcol is None:
            raise ValueError("Feature table needs ImageName or File_Name")
        missing = set(FEATURE_FILES) - set(table.columns)
        if missing:
            raise ValueError(f"Environmental-feature columns missing: {sorted(missing)}")
        table.index = image_ids(table[idcol])
        if table.index.duplicated().any() or set(table.index) != set(order):
            raise ValueError("Feature table must contain every stimulus exactly once")
        aligned = table.loc[order, list(FEATURE_FILES)].apply(pd.to_numeric, errors="raise")
        for name in FEATURE_FILES:
            rdms[name] = scalar_rdm(aligned[name], standardize=True)
    else:
        for name, filename in FEATURE_FILES.items():
            path = Path(rdm_dir) / filename
            rdms[name] = read_labeled_rdm(path, order)
    vectors = np.column_stack([rdm_vector(rdms[name]) for name in FEATURE_FILES])
    if np.any(np.std(vectors, axis=0, ddof=1) <= np.finfo(float).eps):
        raise ValueError("At least one full-set predictor RDM is constant")
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    labels = [x + ".png" for x in order]
    for name, matrix in rdms.items():
        pd.DataFrame(matrix, index=labels, columns=labels).to_csv(
            destination / FEATURE_FILES[name], index_label="ImageName", encoding="utf-8-sig")
    if aligned is not None:
        aligned.to_csv(destination / "features_aligned.csv", index_label="ImageName")
    pd.DataFrame(np.corrcoef(vectors.T), index=FEATURE_FILES, columns=FEATURE_FILES).to_csv(
        destination / "predictor_rdm_pearson.csv", index_label="Predictor")
    pd.DataFrame({"Predictor": list(FEATURE_FILES), "VIF": compute_vif(_zscore_columns(vectors))}).to_csv(
        destination / "predictor_rdm_vif.csv", index=False)


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--feature-csv", type=Path)
    source.add_argument("--rdm-dir", type=Path)
    parser.add_argument("--image-manifest", type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args(argv)
    prepare_feature_rdms(output_dir=args.out_dir, feature_csv=args.feature_csv,
                         rdm_dir=args.rdm_dir, image_manifest=args.image_manifest)


if __name__ == "__main__":
    main()
