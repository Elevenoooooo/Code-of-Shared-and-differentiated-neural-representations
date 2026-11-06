"""Study labels and explicitly scoped Benjamini-Hochberg correction."""

from itertools import combinations
import numpy as np

TASKS = ("Safety", "Vitality", "Beauty", "Imageability")
TASK_PAIRS = tuple(combinations(TASKS, 2))


def canonical_task(value):
    """Return a canonical task label; reject unknown aliases."""
    aliases = {
        "safety": "Safety", "sa": "Safety", "saf": "Safety",
        "vitality": "Vitality", "vi": "Vitality", "vit": "Vitality", "li": "Vitality",
        "beauty": "Beauty", "be": "Beauty", "bea": "Beauty",
        "imageability": "Imageability", "im": "Imageability", "img": "Imageability",
    }
    key = str(value).strip().casefold()
    if key not in aliases:
        raise ValueError(f"Unknown task label: {value!r}")
    return aliases[key]


def fdr_bh(pvalues):
    """Return BH-adjusted P values for one explicitly supplied test family.

    NaNs remain NaN and are not counted as tested hypotheses. The caller must
    choose the correction family; distinct outcome families are not pooled.
    """
    p = np.asarray(pvalues, dtype=float)
    flat = p.ravel()
    valid = np.isfinite(flat)
    if np.any((flat[valid] < 0) | (flat[valid] > 1)):
        raise ValueError("Finite P values must lie in [0, 1].")
    result = np.full(flat.shape, np.nan)
    values = flat[valid]
    if values.size:
        order = np.argsort(values, kind="stable")
        ranked = values[order] * values.size / np.arange(1, values.size + 1)
        adjusted = np.minimum.accumulate(ranked[::-1])[::-1]
        corrected = np.empty_like(adjusted)
        corrected[order] = np.minimum(adjusted, 1)
        result[valid] = corrected
    return result.reshape(p.shape)
