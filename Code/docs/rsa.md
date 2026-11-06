# Representational similarity analysis

This module implements the environmental-feature and perceptual-rating analyses in Figure 5. The four tasks are Safety, Vitality, Beauty, and Imageability. Each participant contributes 99 retained stimulus trials per run, yielding a 99 × 99 neural RDM for each ROI.

## Software and configuration

The MATLAB functions require MATLAB R2022b, Statistics and Machine Learning Toolbox, and SPM12 for NIfTI voxel extraction. The Python module requires NumPy, pandas, and SciPy.

Paths are provided through `study_config`:

| Configuration field | Input or destination |
|---|---|
| `conditions_root` | Task condition files and stimulus metadata |
| `lss_root` | Single-trial LSS beta images |
| `lss_run_manifest` | Manifest associating LSS runs with trial tables |
| `rsa_run_manifest` | Manifest listing extracted ROI pattern files |
| `hcpex_atlas` | Bilateral merged HCPex17 label NIfTI |
| `hcpex_labels` | Atlas ROI identifiers and labels |
| `spm_dir` | SPM12 installation |
| `image_id_csv` | Stimulus manifest containing `File_Name` and `Img_ID` |
| `urban_rdm_dir` | Ten environmental-feature RDM CSV files |
| `derivatives_root` | Extracted patterns and trial tables |
| `output_root` | RSA coefficients, statistics, and quality-control tables |

Stimulus filenames use `image1.png` through `image400.png`. `File_Name` contains this standard stimulus label, while `Img_ID` is the experiment index. Image labels identify stimuli when aligning different input tables.

## Environmental-feature RDMs

The simultaneous model contains ten predictors:

| Feature | RDM filename |
|---|---|
| Colorfulness | `RDM_colorfulness.csv` |
| Visual Complexity | `RDM_visual_complexity.csv` |
| Spatial Layering | `RDM_spatial_layering.csv` |
| Green View Index | `RDM_green_view_index.csv` |
| Sky Visibility | `RDM_sky_visibility.csv` |
| Enclosure | `RDM_enclosure.csv` |
| Crowd Attraction | `RDM_crowd_attraction.csv` |
| Vehicle Interference | `RDM_vehicle_interference.csv` |
| Street Facility Convenience | `RDM_street_facility_convenience.csv` |
| Urbanization Intensity Index | `RDM_urbanization_intensity_index.csv` |

For feature k, scalar measurements are standardized across all 400 stimuli using their mean and sample standard deviation. Feature dissimilarity between images i and j is the absolute difference between their standardized values:

\[
d_{ij}^{(k)} = |z_{ik}-z_{jk}|.
\]

The feature table contains one row per stimulus, an `ImageName` or `File_Name` column, and ten numeric columns matching the names above. Input validation rejects missing or constant measurements. From the code root:

```powershell
python -m urban_perception.rsa --feature-csv <features.csv> --out-dir <derivatives/urban_rdms>
```

Labeled 400 × 400 RDMs can also be supplied directly:

```powershell
python -m urban_perception.rsa --rdm-dir <feature_rdms> --out-dir <derivatives/urban_rdms>
```

Each CSV contains a first image-label column followed by 400 image-label columns. Row and column labels are aligned independently. Matrices must contain all 400 unique stimuli, finite nonnegative distances, a zero diagonal, and symmetric entries within an absolute tolerance of 1e-8.

The output order is `image1.png` through `image400.png`. The optional argument `--image-manifest config/image_id.csv` uses experiment-index order instead. The output includes the ten aligned RDMs, Pearson correlations among their vectorized distances, and variance inflation factors (VIF). When the input is a scalar feature table, `features_aligned.csv` contains the aligned measurements.

## Single-trial patterns and manifests

The MATLAB workflow is:

```matlab
setup_study_paths();
cfg = study_config();
build_rsa_manifest(cfg);
extract_rsa_patterns(cfg);
run_rsa(cfg);
```

`build_rsa_manifest` reads task-condition `StimMeta_RT` tables containing `BetaIndex`, `OriginalTrial`, `ImageName`, and `KeyRaw`. Each table represents 99 retained stimuli. These fields associate each beta image with its stimulus and rating. LSS beta filenames follow `beta_stim_trial_001.nii` through `beta_stim_trial_099.nii`.

The generated LSS manifest contains:

| Column | Description |
|---|---|
| `Group` | Latin-square group, `Group1`–`Group4` |
| `Subject` | Participant identifier |
| `Run` | Run identifier, `Run1`–`Run4` |
| `Task` | Perceptual task |
| `BetaDir` | Folder of single-trial LSS beta images |
| `TrialTableFile` | CSV with `BetaFile`, `ImageName`, and `Rating` |

`extract_rsa_patterns` samples atlas labels in each run's beta grid using nearest-neighbour interpolation and extracts voxel values independently for each ROI. Trial tables define beta order, and duplicate beta files or stimulus identities are rejected. Optional `RT` values are retained as metadata.

The resulting MAT files contain `PatternCell`, `roi_labels`, `roi_ids`, `pathRaw`, and `ratingVec`. Each `PatternCell{roi}` is a 99 × voxel matrix. `pathRaw` and `ratingVec` follow the same trial order. The RSA manifest lists `Group`, `Subject`, `Run`, `Task`, `RunPatternFile`, and `RatingTableFile`. Relative paths are resolved against `derivatives_root`.

Pattern files can identify stimuli using `pathRaw` or `ImageNames`. When `canonIdx` is supplied, its experiment indices are joined to the `Img_ID` column in `image_id_csv`; the corresponding `File_Name` supplies each stimulus label. A separate rating table contains `ImageName` and `Rating` if the pattern file does not include `ratingVec`.

`make_rsa_trial_table` provides a behavior-table interface using the selected rating column and the number of formal trials preceding the retained beta sequence. Its output records beta filenames, standard stimulus labels, ratings, and behavior row/trial indices.

## Neural dissimilarity and simultaneous RSA–GLM

Neural dissimilarity is one minus the Pearson correlation between two trials' voxel patterns within an ROI:

\[
d_{ij}^{\mathrm{neural}}=1-r(V_i,V_j).
\]

Finite, nonconstant voxels contribute to the pattern correlation. Each model uses 99 retained trials and 4,851 unique image pairs. Amygdala and Entorhinal cortex are excluded, retaining 15 ROIs. The environmental-feature GLM uses VVC, PPC (the merged SPL/IPL ROI), and vmPFC.

For each participant, task, and selected ROI, the 400-image feature RDMs are subset to the 99 observed stimulus labels in neural-pattern order. Strict lower-triangular vectors are extracted from all matrices using the same pair ordering. The neural response vector and each feature vector are separately standardized, and an intercept and all ten predictors enter one OLS model:

\[
\mathbf y=\beta_0+\sum_{k=1}^{10}\beta_k\mathbf X_k+\boldsymbol\varepsilon.
\]

Coefficients describe the unique correspondence of each feature with neural representational geometry while controlling for the other features. Outputs include standardized coefficients, VIF, model R², adjusted R², and each predictor's unique contribution to explained variance. A constant within-run predictor is recorded in the model-fit table and omitted from that fit.

Run estimates are averaged within participant, task, ROI, and predictor. Two-sided one-sample t tests evaluate participant coefficients against zero. Benjamini–Hochberg FDR correction is applied across the three selected ROIs, four tasks, and ten features: 120 tests when all model cells are available. The primary adjusted value is `Q_FDR`, with significance defined as q < 0.05.

## Perceptual-rating RSA

Each participant's task-specific rating RDM contains absolute response differences, d(i,j) = |r(i) − r(j)|. Spearman correlations between rating and neural RDM vectors are calculated in all 15 retained ROIs. Pairs involving missing responses are omitted. Correlations are clipped to the interval [−0.999999, 0.999999] before Fisher transformation.

Group tests use participant Fisher-z values in two-sided one-sample t tests against zero. Benjamini–Hochberg FDR correction is applied within each ROI across the four tasks. `MeanSpearmanRho` is the arithmetic mean of participant Spearman correlations, and `Q_FDR` contains the adjusted group-test p value.

## Outputs

Results are written to `<output_root>/rsa`:

| File | Contents |
|---|---|
| `subject_run_feature_beta.csv` | Standardized feature coefficients, VIF, and unique explained variance |
| `subject_run_rating_rsa.csv` | Spearman correlations, Fisher z values, and pair counts |
| `group_feature_beta.csv` | Group means, SEM, t tests, and FDR-adjusted values |
| `group_rating_rsa.csv` | Group rating-RSA estimates and FDR-adjusted values |
| `model_fit.csv` | Model fit and constant-predictor records |
| `neural_rdm_qc.csv` | Retained trial/voxel counts and extraction status |
| `run_manifest_inclusion.csv` | Participant/run inclusion records |
| `predictor_rdm_pearson.csv` | Correlations among feature-RDM vectors |
| `predictor_rdm_vif.csv` | Full-stimulus-set VIF values |
| `rsa_results.mat` | Result tables and analysis configuration |

The default analysis requires complete cases across all four tasks and enforces 99 retained trials per neural RDM. Quality-control thresholds are 20 valid trials, 10 finite nonconstant voxels, and 1,000 finite pairs; the 99-trial requirement applies in addition to these minimum thresholds. Feature group tests require at least eight participants, and rating group tests require at least two.
