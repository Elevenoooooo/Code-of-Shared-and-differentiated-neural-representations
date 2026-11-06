# Data formats and configuration

## Companion data layout

```text
Data/
  Images/                  image1.png ... image400.png
  Behavior/
    Group1/Run1/ ...        behavioural trial CSVs
    Group2/Run1/ ...
    Group3/Run1/ ...
    Group4/Run1/ ...
  ETHeatmap/               task-specific fixation-map visualizations
  fMRI/                    group statistical maps
```

Additional derivatives include preprocessed time series, motion parameters, event MAT files, single-trial beta images, ROI voxel patterns, numeric gaze maps and model RDMs. Paths are specified through `study_config.m` or the Python command-line arguments.

## Participants, tasks and stimuli

`config/participants.csv` has `Group,Subject,Include` columns. Subject is a zero-padded participant identifier. Group and Subject jointly identify a participant. The MATLAB list `config/subjects.txt` organizes these identifiers by group and run.

The task schedule is:

| Group | Run 1 | Run 2 | Run 3 | Run 4 |
|---|---|---|---|---|
| 1 | Safety | Vitality | Beauty | Imageability |
| 2 | Vitality | Safety | Imageability | Beauty |
| 3 | Beauty | Imageability | Safety | Vitality |
| 4 | Imageability | Beauty | Vitality | Safety |

`config/image_id.csv` contains `File_Name,Img_ID`. File_Name identifies the stimulus, and Img_ID records its experimental index. Matrices are aligned by these identifiers. Behavioural trial files contain `Path_raw`, `thisTrialN_raw`, response keys and response times. `Path_raw` identifies the stimulus presented on that trial.

## Functional imaging

The preprocessed input hierarchy follows `GroupN/RunN/<functional-folder>/UP_<subject>/`. A run contains a four-dimensional NIfTI and matching realignment parameters. Motion files follow `rp_aUP_<subject>_BOLD_PA_TASK<run>_*.txt`. Configured NIfTI basenames include `wFiltered_4DVolume.nii`, `w4DVolume.nii` and `w4D.nii`.

Event files are stored in `conditions/conds_glm_GroupN_RunN_RT/`. Conditions are Stim, Rating and Fixation. Each retained event has an onset and duration; the conventional GLM additionally uses trial-order and response pmods. StimMeta CSVs contain BetaIndex, OriginalTrial, ImageName, Img_ID and response metadata. There are 99 retained trials per run.

LSS stimulus beta files use `beta_stim_trial_###.nii` under `lss/GroupN/RunN/UP_<subject>/LSS_StimBeta_RunN/`.

Group outputs include unthresholded statistics, p/q maps and q<0.05 significance masks. The contrast prefix is `Stim_gt_Fixation`. Statistical images and atlas labels must use matching spatial coordinates for voxel summaries.

## Gaze tables

The trial-QC table contains `Group,Subject,Run,Trial,Task,ImageName,CalibrationOK,GazeValidRatio,ValidRecordingDuration_ms`.

The fixation table contains `Group,Subject,Run,Trial,X,Y,Start_ms,End_ms`. Coordinates are pixels and event times are milliseconds relative to stimulus onset. Analysis uses the four-second stimulus window. Numeric group-map indices contain `Task,ImageName,HeatmapFile`; supported probability-map formats are detailed in the gaze methods.

## Feature and neural RDMs

The feature table has one row per stimulus and the ten environmental predictors specified in the RSA methods. A labelled feature RDM has a stimulus-identifier column and matching stimulus columns. The full stimulus model has 400 rows/columns; each run selects the corresponding 99 x 99 submatrix.

The LSS manifest contains `Group,Subject,Run,Task,BetaDir,TrialTableFile`. The RSA manifest contains `Group,Subject,Run,Task,RunPatternFile`. Relative manifest paths are resolved under cfg.derivatives_root. Pattern MATs contain the trial-by-voxel patterns, ROI labels, ordered image identifiers and task-specific ratings. See [RSA](rsa.md) for field definitions.

## Configuration

`study_config.m` sets data, derivative, software, atlas and output locations. Acquisition-specific slice timing and input sessions are supplied in the DPABI Cfg MAT. Core parameters are also listed in `config/study_parameters.json`.
