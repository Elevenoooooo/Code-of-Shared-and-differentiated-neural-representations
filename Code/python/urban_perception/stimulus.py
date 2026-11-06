"""Compare feature distributions of a stimulus subset and its parent pool."""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2, ks_2samp, shapiro

from .behavior import normalize_image_name


def distribution_diagnostics(parent: pd.DataFrame, selected: pd.DataFrame,
                              features: list[str], image_column: str = "File_Name"
                              ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Summarize univariate distributions, PCA coverage and Mahalanobis distances.

    Selected images must be a true subset and carry identical feature values.
    Missing/nonfinite features fail rather than silently reducing the sample.
    KS P values are descriptive: selected and parent sets are not independent.
    """
    for label, frame in (("parent", parent), ("selected", selected)):
        if not set(features + [image_column]).issubset(frame):
            raise ValueError(f"Missing {label} feature or image columns")
        if frame[image_column].isna().any() or frame[image_column].duplicated().any():
            raise ValueError(f"Missing/duplicate {label} image identifiers")
        if len(frame) < 3:
            raise ValueError(f"At least three {label} images required")
    parent, selected = parent.copy(), selected.copy()
    for frame in (parent, selected):
        frame[image_column] = frame[image_column].map(normalize_image_name)
        if frame[image_column].duplicated().any():
            raise ValueError("Duplicate image identifiers")
    p = parent.set_index(image_column)[features].apply(pd.to_numeric, errors="raise")
    s = selected.set_index(image_column)[features].apply(pd.to_numeric, errors="raise")
    if not np.isfinite(p.to_numpy()).all() or not np.isfinite(s.to_numpy()).all():
        raise ValueError("Feature values must be finite")
    if not s.index.isin(p.index).all():
        raise ValueError("Selected IDs must be a subset of parent IDs")
    if not np.allclose(p.loc[s.index].to_numpy(), s.to_numpy(), rtol=1e-7, atol=1e-10):
        raise ValueError("Parent and selected rows disagree on their feature values")
    rows = []
    for name in features:
        x, y = p[name].to_numpy(), s[name].to_numpy()
        pooled_sd = np.sqrt(((len(x)-1)*x.var(ddof=1) + (len(y)-1)*y.var(ddof=1)) / (len(x)+len(y)-2))
        low, high = np.percentile(x, [10, 90])
        ks = ks_2samp(x, y)
        rows.append(dict(Feature=name, N_Parent=len(x), N_Selected=len(y),
                         Mean_Parent=x.mean(), Mean_Selected=y.mean(),
                         SD_Parent=x.std(ddof=1), SD_Selected=y.std(ddof=1),
                         Cohen_d=(y.mean()-x.mean())/pooled_sd if pooled_sd > 0 else np.nan,
                         KS_Statistic=ks.statistic, KS_P_descriptive=ks.pvalue,
                         Shapiro_P_Selected=shapiro(y).pvalue if len(y) <= 5000 and y.std() > 0 else np.nan,
                         Parent_P10=low, Parent_P90=high, Selected_Min=y.min(), Selected_Max=y.max(),
                         Covers_Parent_P10_P90=bool(y.min() <= low and y.max() >= high),
                         Selected_Outside_P10_P90=float(((y < low) | (y > high)).mean())))
    # Fit standardization and principal axes to the parent feature distribution.
    mu, sd = p.mean().to_numpy(), p.std(ddof=0).to_numpy()
    sd[sd == 0] = 1
    x, y = (p.to_numpy()-mu)/sd, (s.to_numpy()-mu)/sd
    _, singular, vt = np.linalg.svd(x, full_matrices=False)
    if len(singular) < 2 or np.square(singular).sum() == 0:
        raise ValueError("At least two informative dimensions required for PCA diagnostics")
    explained = singular**2 / np.square(singular).sum()
    coordinates = []
    for label, ids, values in (("Parent", p.index, x), ("Selected", s.index, y)):
        scores = values @ vt[:2].T
        for image, score in zip(ids, scores):
            coordinates.append({image_column: image, "Set": label, "PC1": score[0], "PC2": score[1],
                                "PC1_ExplainedVarianceRatio": explained[0],
                                "PC2_ExplainedVarianceRatio": explained[1]})
    covariance = np.atleast_2d(np.cov(x, rowvar=False)) + np.eye(len(features))*1e-8
    centered = y - x.mean(axis=0)
    distances = np.sqrt(np.maximum(0, np.einsum("ij,jk,ik->i", centered, np.linalg.inv(covariance), centered)))
    threshold = np.sqrt(chi2.ppf(.99, len(features)))
    outliers = pd.DataFrame({image_column: s.index, "Mahalanobis": distances,
                             "ChiSquareApprox_99PercentThreshold": threshold,
                             "BeyondApproxThreshold": distances > threshold})
    return pd.DataFrame(rows), pd.DataFrame(coordinates), outliers


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    compare = commands.add_parser("diagnostics")
    compare.add_argument("--parent-features", type=Path, required=True)
    compare.add_argument("--selected-features", type=Path, required=True)
    compare.add_argument("--image-column", default="File_Name")
    compare.add_argument("--features", nargs="+", required=True)
    compare.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    result, pca, outliers = distribution_diagnostics(pd.read_csv(args.parent_features),
        pd.read_csv(args.selected_features), args.features, args.image_column)
    result.to_csv(args.output_dir / "stimulus_distribution_diagnostics.csv", index=False)
    pca.to_csv(args.output_dir / "stimulus_pca_coordinates.csv", index=False)
    outliers.to_csv(args.output_dir / "stimulus_mahalanobis_diagnostics.csv", index=False)


if __name__ == "__main__":
    main()
