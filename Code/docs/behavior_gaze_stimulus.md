# Behaviour, gaze and stimulus analyses

The Python modules provide image-level behavioural summaries, fixation-map similarity analyses and stimulus-feature summaries. Stimuli use the identifiers `image1.png` to `image400.png`. The four perceptual tasks are Safety, Vitality, Beauty and Imageability.

Commands below are relative to the `Code` directory. Each analysis writes its tables to the specified output directory.

## Behavioural ratings

### Inputs

The behavioural directory contains one CSV per participant and run, named `GroupN_NNN_runN_<task>_trials.csv`. Each participant is identified by the combination of `Group` and `Subject`.

| Field | Definition |
| --- | --- |
| `Path_raw` | Path to the image presented on the trial; its basename provides the image identifier |
| `thisTrialN_raw` | Zero-based formal trial number |
| `key_resp_sc.keys_raw` | Binary response key: 1 for the negative response and 2 for the positive response |
| `key_resp_sc.rt_raw` | Response time in seconds |

Task abbreviations in the CSV filenames are `sa` for Safety, `li` for Vitality, `bea` for Beauty and `img` for Imageability. The filename task is checked against the run schedule:

| Group | Run 1 | Run 2 | Run 3 | Run 4 |
| --- | --- | --- | --- | --- |
| Group1 | Safety | Vitality | Beauty | Imageability |
| Group2 | Vitality | Safety | Imageability | Beauty |
| Group3 | Beauty | Imageability | Safety | Vitality |
| Group4 | Imageability | Beauty | Vitality | Safety |

The inclusion table is `config/participants.csv`, with columns `Group,Subject,Include`. `Include` accepts `true/false`, `yes/no` or `1/0`. An optional `Run` column supports run-specific inclusion; a blank or omitted Run applies to all four runs.

### Method

Formal trials are selected using numeric `thisTrialN_raw` values. Response keys 1 and 2 are coded as -1 and +1, respectively. Other keys are treated as missing. Ratings are averaged within participant, task and image before averaging across participants with equal weight.

For each task and image, the outputs include the net rating mean, standard deviation and number of participants. `PositiveRate=(NetRatingMean+1)/2`; `NetScore_100=100*NetRatingMean` ranges from -100 to +100.

Pearson correlations are calculated across corresponding image-level rating profiles for the six task pairs. Benjamini-Hochberg correction is applied across these six comparisons. Each comparison reports its number of complete image pairs.

### Command and outputs

```console
python -m urban_perception.behavior binary --behavior-dir ../Data/Behavior --participants config/participants.csv --output-dir outputs/behavior
```

| Output | Contents |
| --- | --- |
| `binary_rating_trials.csv` | Trial identifiers, image identifiers, raw keys, signed ratings and response times |
| `behavior_qc.csv` | Inclusion status, formal-trial counts, valid-response counts and repeated-image counts by run |
| `image_level_rating_profiles.csv` | Task-by-image rating summaries |
| `rating_profiles_wide.csv` | Image-by-task positive-response rates |
| `rating_correlations.csv` | Six Pearson correlations, P values, BH-adjusted P values and image counts |
| `analysis_metadata.json` | Analysis parameters and participant, run and image counts |

### Online image-level ratings

The `online` command accepts one row per image with a unique `imageN.png` identifier and the columns `safety,liveliness,beauty,imageability`. Canonical task labels are also accepted. Ratings in this table are image-level means.

```console
python -m urban_perception.behavior online --ratings inputs/online_image_mean_ratings.csv --image-column File_Name --output-dir outputs/online_ratings
```

Outputs are `rating_profiles_wide.csv`, `rating_correlations.csv` and `analysis_metadata.json`.

## Fixation maps and cross-task gaze similarity

### Trial and fixation inputs

The trial table contains one row per trial:

```csv
Group,Subject,Run,Trial,Task,ImageName,CalibrationOK,GazeValidRatio,ValidRecordingDuration_ms
Group1,001,1,1,Safety,image1.png,true,0.95,4000
```

`CalibrationOK` indicates calibration success. `GazeValidRatio` is the proportion of valid gaze samples in the 4-s image-viewing window. `ValidRecordingDuration_ms` gives the valid recording duration for that window. Trials are included when calibration succeeds, the valid-sample ratio is at least 0.70 and recording duration is at least 3,900 ms.

The fixation table contains EyeLink-defined fixation events:

```csv
Group,Subject,Run,Trial,X,Y,Start_ms,End_ms
Group1,001,1,1,960.0,540.0,100.0,350.0
```

`X,Y` are zero-based pixels on a 1,920 by 1,080 display. `Start_ms,End_ms` are event times relative to image onset. The columns `Group,Subject,Run,Trial` link fixation events to their trials.

### Probability-map construction

Fixation events are clipped to the interval from 0 to 4,000 ms and weighted by their duration within that interval. Coordinates are accumulated in 4-pixel bins, producing a 270 by 480 map. A Gaussian kernel with standard deviation 6 bins and radius `ceil(3*sigma)` provides spatial smoothing; convolution uses zero padding. Each trial map is normalized to unit sum.

Maps are averaged within participant before averaging participants with equal weight for each task and image. This yields one group map per task and image.

```console
python -m urban_perception.gaze from-fixations --trials inputs/gaze_trials.csv --fixations inputs/fixations.csv --participants config/participants.csv --output-dir outputs/gaze_numeric
```

Outputs are the numeric `.npy` maps in `numeric_maps/`, `group_map_index.csv`, `gaze_trial_qc.csv` and `analysis_metadata.json`.

### Numeric map index

The similarity analysis accepts a table with one row per task and image:

```csv
Task,ImageName,HeatmapFile
Safety,image1.png,numeric_maps/Safety_image1.npy
```

Map files are nonnegative two-dimensional arrays with positive mass. Supported formats are `.npy`, `.npz` containing the array `H`, and MATLAB `.mat` containing `H`. All maps share the same dimensions. Relative file paths are resolved against the index table's directory. The default image set is `image1.png` through `image400.png` for each of the four tasks.

### Similarity statistics

Pearson correlation pools corresponding map locations across images. Structural similarity (SSIM) is computed separately for each image and then averaged across images. SSIM uses the paired maps' shared data range and an odd window size of up to 7 pixels.

Whole images are resampled with replacement 5,000 times, using seed 20260717. The analysis reports percentile 95% confidence intervals and the one-sided empirical tail probability `(number of bootstrap estimates <= 0 + 1)/(number of bootstrap samples + 1)` for positive similarity. Benjamini-Hochberg correction is applied across the six task pairs separately for Pearson correlation and SSIM.

```console
python -m urban_perception.gaze similarity --index outputs/gaze_numeric/group_map_index.csv --n-bootstrap 5000 --seed 20260717 --output-dir outputs/gaze_similarity
```

| Output | Contents |
| --- | --- |
| `gaze_task_similarity.csv` | Pearson and SSIM estimates, confidence intervals, bootstrap tail probabilities and BH-adjusted P values for the six task pairs |
| `image_level_ssim.csv` | Image-specific SSIM values for every task pair |
| `analysis_metadata.json` | Input index, image count, bootstrap count and random seed |

## Stimulus-feature summaries

### Low-level visual features

`urban_perception.features` calculates brightness, colourfulness, GLCM texture contrast, GLCM entropy and grayscale-image entropy. Brightness is the mean grayscale intensity divided by 255. Colourfulness is calculated after resizing to a width of 250 pixels while preserving aspect ratio. GLCMs use a distance of one pixel and orientations of 0, 45, 90 and 135 degrees. Contrast is averaged across orientations; entropy is summed across the four orientation-specific probability matrices.

```console
python -m urban_perception.features --images-dir ../Data/Images --output outputs/stimulus/rgb_features.csv
```

The output has one row per image, with columns `File_Name,brightness,colorfulness,texture_contrast_glcm,texture_entropy_glcm,visual_entropy_image`. Feature extraction uses `opencv-python` in addition to the core Python dependencies.

The optional `--depth-dir` argument accepts quantitative single-channel `.npy` depth arrays with the same image stem. Depth summaries are the standard deviation, mean Sobel-gradient magnitude and absolute difference between the mean values of the left and right halves. Their output columns are `depth_std,depth_gradient_mean,depth_symmetry_diff`.

### Feature-distribution comparisons

`urban_perception.stimulus` summarizes the selected stimulus set relative to a parent image pool. Both tables contain unique image identifiers and the numeric feature columns specified with `--features`. Selected rows share their identifiers and feature values with the corresponding parent rows.

```console
python -m urban_perception.stimulus diagnostics --parent-features inputs/parent_features.csv --selected-features inputs/selected_features.csv --image-column File_Name --features brightness colorfulness texture_contrast_glcm visual_entropy_image depth_std --output-dir outputs/stimulus_distributions
```

Outputs include descriptive means and standard deviations, Cohen's d, descriptive KS statistics, Shapiro-Wilk P values and parent P10-P90 coverage. PCA uses standardization and principal axes fitted to the parent pool. Mahalanobis distances use the parent covariance with a diagonal regularization term of 1e-8; the reported reference threshold is the 99th percentile of the chi-square approximation.

The three output tables are `stimulus_distribution_diagnostics.csv`, `stimulus_pca_coordinates.csv` and `stimulus_mahalanobis_diagnostics.csv`.
