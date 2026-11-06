"""Low-level visual features of street-view stimuli and quantitative depth maps."""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image
from scipy.ndimage import sobel
from skimage.feature import graycomatrix, graycoprops

from .behavior import normalize_image_name


def rgb_features(rgb: np.ndarray) -> dict[str, float]:
    """Calculate grayscale, colourfulness and texture features from an RGB array."""
    import cv2
    if rgb.ndim != 3 or rgb.shape[2] != 3 or rgb.dtype != np.uint8:
        raise ValueError("Expected an 8-bit RGB array")
    gray = cv2.cvtColor(rgb, cv2.COLOR_RGB2GRAY)
    # Calculate colourfulness at a width of 250 pixels, preserving aspect ratio.
    small = cv2.resize(rgb, (250, int(rgb.shape[0] * 250 / rgb.shape[1])), interpolation=cv2.INTER_AREA)
    r, g, b = small.astype(float).transpose(2, 0, 1)
    rg, yb = np.abs(r-g), np.abs((r+g)/2-b)
    colour = np.hypot(rg.std(), yb.std()) + 0.3*np.hypot(rg.mean(), yb.mean())
    glcm = graycomatrix(gray, distances=[1], angles=[0, np.pi/4, np.pi/2, 3*np.pi/4],
                       levels=256, symmetric=True, normed=True)
    hist = np.bincount(gray.ravel(), minlength=256).astype(float)
    probabilities = hist / hist.sum()
    # Sum entropy over four separately normalized orientation-specific GLCMs.
    return {"brightness": float((gray.astype(np.float32)/255).mean()),
            "colorfulness": float(colour),
            "texture_contrast_glcm": float(graycoprops(glcm, "contrast").mean()),
            "texture_entropy_glcm": float(-np.sum(glcm*np.log2(glcm+(glcm == 0)))),
            "visual_entropy_image": float(-np.sum(probabilities*np.log2(probabilities+1e-12)))}


def depth_features(depth: np.ndarray) -> dict[str, float]:
    """Calculate depth variation, Sobel gradient and left-right mean difference."""
    array = np.asarray(depth, dtype=float)
    if array.ndim != 2 or array.shape[1] < 2 or not np.isfinite(array).all():
        raise ValueError("Expected finite 2-D quantitative depth with width >=2")
    gradient = np.hypot(sobel(array, axis=1), sobel(array, axis=0))
    middle = array.shape[1] // 2
    return {"depth_std": float(array.std()), "depth_gradient_mean": float(gradient.mean()),
            "depth_symmetry_diff": float(abs(array[:, :middle].mean()-array[:, middle:].mean()))}


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--images-dir", type=Path, required=True)
    parser.add_argument("--depth-dir", type=Path, help="Optional same-stem .npy quantitative depth arrays")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    files = sorted(p for p in args.images_dir.iterdir() if p.suffix.lower() == ".png")
    if not files:
        raise FileNotFoundError("No RGB image files found")
    results = []
    for path in files:
        with Image.open(path) as image:
            rgb = np.array(image.convert("RGB"))
        row = {"File_Name": normalize_image_name(path.name), **rgb_features(rgb)}
        if args.depth_dir:
            depth_path = args.depth_dir / (path.stem + ".npy")
            if not depth_path.exists():
                raise FileNotFoundError(f"Requested quantitative depth missing: {depth_path}")
            row.update(depth_features(np.load(depth_path, allow_pickle=False)))
        results.append(row)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(results).to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
