# Activation and functional connectivity

This module implements the atlas, event, activation and beta-series connectivity analyses for *Shared and differentiated neural representations of urban scenes in the human brain*. All functions accept the common structure returned by `study_config`.

## Analysis functions

| Function | Purpose | Principal output |
|---|---|---|
| `build_merged_hcpex_atlas(cfg)` | Merge HCPex parcels into 17 bilateral ROIs | Label atlas, parcel lookup and voxel counts |
| `update_atlas_voxel_counts(cfg)` | Add ROI denominators to the derivative lookup | Updated lookup path |
| `run_build_conditions(cfg)` | Construct trial events and stimulus parametric modulators | Condition MATs, trial metadata and valid-response subject list |
| `run_activation_first_level(cfg)` | Fit each participant/run activation GLM | SPM model and Stim > Fixation contrast |
| `run_activation_group_model(cfg)` | Test positive activation and differences between tasks | One-sample T maps, omnibus F maps, signed pairwise T maps and FDR products |
| `run_activation_conjunction(cfg)` | Compute common activation across all four tasks | Intersection mask, minimum-T map and ROI summaries |
| `split_activation_overlap_masks(cfg)` | Export task-overlap masks for display | Individual bit-code masks and legend |
| `run_lss_first_level(cfg)` | Estimate each retained stimulus trial separately | 99 target-trial beta images per run |
| `run_beta_series_connectivity(cfg)` | Extract ROI beta series and calculate connectivity | Individual correlation/z matrices and complete-case edge tables |
| `fit_connectivity_task_lme(cfg, EdgeLong)` | Test the omnibus task effect on each connection | Edge-level F, degrees of freedom, P and FDR-adjusted Q |

The `EdgeLong` argument is optional: without it, `fit_connectivity_task_lme` reads `BSC_EdgeLong.csv` from `cfg.connectivity_root`.

## Regions of interest

HCPex parcels are combined according to the workbook's `Parcel_Map` and `ROI_Definition` sheets. Each source parcel belongs to exactly one merged bilateral ROI. The atlas outputs preserve integer labels and list the constituent parcels and voxel counts for each region.

Subsequent analyses retain 15 regions by excluding the amygdala and entorhinal cortex. ROI labels are sampled in functional-image space by nearest-neighbour interpolation. The activation group model uses the first available participant contrast as a geometry reference and writes `Atlas/HCPex_Merged17_functional.nii`. The ROI union and voxel denominators are calculated in that same functional grid. Participant contrasts already on this grid retain their spatial sampling. Beta-series extraction likewise samples atlas labels in each beta-image grid. Image geometry includes both voxel dimensions and the voxel-to-world affine.

## Events and first-level activation

Each run contains 99 retained trials following removal of the first six-second trial together with the initial six volumes at TR = 1 s. The condition builder reads the behavioural table, response time and response key for each observed trial. The event names are:

- **Stim:** the four-second streetscape presentation.
- **Rating:** the response interval, beginning after presentation and lasting until the recorded response.
- **Fixation:** the response-following fixation period, occupying the remainder of the two-second response/fixation interval.

Valid response times lie strictly between zero and two seconds. Trial order is standardized within each run and entered as a first-order Stim parametric modulator. The response modulator codes key 1 as -1, key 2 as +1 and missing/invalid responses as zero. Sequential orthogonalization of the two stimulus modulators is disabled.

The activation GLM convolves events with the canonical haemodynamic response function. Nuisance regressors contain Friston-24 expanded realignment terms and one indicator per high-motion volume. Framewise displacement uses a threshold of 0.5 mm and a 50-mm radius to convert rotations to displacement. Other defaults are TR = 1 s, microtime resolution/reference = 60/30, high-pass cutoff = 128 s, no HRF derivatives, AR(1) serial correlations and a run intercept.

The primary contrast is **Stim > Fixation**, with weights +1 for Stim, -1 for Fixation and zero for all other regressors, including parametric modulators. Its output stem is `Stim_gt_Fixation`.

The event MAT files contain `names`, `onsets`, `durations`, `pmod` and `orth`. A corresponding `StimMeta_RT_GroupN_subject_RunN.csv` records `BetaIndex`, `OriginalTrial`, `ImageName`, `Img_ID`, `StimOnset`, `OrderParam`, `KeyRaw` and `RespBin`. These columns provide the trial identities used by the RSA manifest builder.

## Group activation and task differentiation

For each task, positive activation is tested using a pooled one-sample T model containing one Stim > Fixation contrast per participant and an intercept. The same complete-case cohort, with all four tasks, contributes to the repeated-measures analysis. The four positive activation maps receive voxelwise Benjamini-Hochberg FDR correction at q < 0.05 across the atlas ROI union.

Task differentiation uses a repeated-measures flexible-factorial design with Subject and a 16-level GroupTaskCell factor. The latter represents all combinations of the four Latin-square groups and four tasks. Subject intercepts absorb participant baselines; task marginal contrasts average equally over groups. This parameterization represents the group-by-task assignment while maintaining estimable task comparisons.

An omnibus F test evaluates equality of the four task marginal responses. Six signed A-minus-B T contrasts describe pairwise differences:

1. Safety - Vitality
2. Safety - Beauty
3. Safety - Imageability
4. Vitality - Beauty
5. Vitality - Imageability
6. Beauty - Imageability

Pairwise probabilities are two-sided. For omnibus and pairwise differentiation maps, voxelwise BH-FDR correction is applied within each ROI. The output retains unthresholded statistics alongside p maps, q maps, significant-voxel masks and ROI summaries. ROI summaries include the mean statistic and fraction of significant voxels. Signed pairwise statistics and masks provide inputs to significance-weighted skewness summaries.

Common activation is the strict intersection of the four positive FDR masks within the final 15 ROIs. The visualization statistic is the minimum unthresholded T across tasks inside this intersection. Overlap bit codes are Safety = 1, Vitality = 2, Beauty = 4 and Imageability = 8; code 15 denotes all four tasks.

## Single-trial beta estimation and connectivity

Each LSS model includes the target Stim trial, all remaining Stim trials, Rating, Fixation, Friston-24 terms, high-motion indicators and an intercept. The target coefficient supplies one beta image. Thus, each retained run contributes 99 target-trial beta images to connectivity and voxelwise RSA.

For connectivity, finite voxel betas are averaged independently within each of the 15 ROIs for each trial, yielding a 99-by-15 matrix. Pearson correlation across trials gives a 15-by-15 connectivity matrix. The 105 unique off-diagonal values are Fisher-z transformed. Correlations are clipped to +/-0.999999 before `atanh`; stored diagonals are zero. Constant or nonfinite ROI series yield unavailable edges.

Connectivity inference uses the model:

```text
Z ~ Task + Group + (1 | Subject)
```

Task and Latin-square Group are fixed effects; Subject has a random intercept. Models use restricted maximum likelihood and effect coding. Omnibus task F tests use residual denominator degrees of freedom. BH-FDR correction is applied across 105 connections at q < 0.05. Complete cases contain all four tasks; participant keys combine group and subject identifiers.

## Configuration

| Field | Input or destination |
|---|---|
| `data_root` | Data root containing `Behavior/GroupN/RunN/` |
| `derivatives_root`, `output_root` | Derived-data and summary destinations |
| `subject_list` | Subject/run inclusion TXT using `GroupN` and `RunN = {'001', ...};` blocks |
| `image_id_csv` | Stimulus table containing `File_Name` and `Img_ID` |
| `n_trials` | Retained trial count; 99 |
| `spm_dir` | SPM12 installation directory |
| `preprocessed_root` | Preprocessed functional volumes and `RealignParameter` folders |
| `conditions_root` | Event MATs, trial metadata and generated subject list |
| `first_level_root` | Participant/run activation models |
| `activation_root` | Group activation models, maps and summaries |
| `lss_root` | LSS models and collected target beta images |
| `connectivity_root` | Beta-series matrices, edge tables and connectivity statistics |
| `hcpex_source_atlas` | HCPex parcel-label NIfTI used for atlas construction |
| `atlas_crosswalk` | Workbook defining the bilateral merged ROIs |
| `merged_atlas_dir` | Merged atlas, lookup and voxel-count destination |
| `hcpex_atlas`, `hcpex_labels` | Merged atlas NIfTI and parcel lookup used by analysis functions |
| `plan` (optional) | Selected `{groupName, {runNames}; ...}` combinations |

Stage-specific options are accepted through `condition_options`, `activation_first_level_options`, `activation_group_options`, `conjunction_options`, `lss_options`, `atlas_options`, `connectivity_options` and `connectivity_lme_options`. The fields correspond to the documented parameters defined near each function's entry point; unknown keys are rejected.

The functional input layout is `GroupN/RunN/<functional-folder>/UP_subject/`; realignment files reside under the corresponding `RealignParameter/UP_subject/` directory. The functional filename priority is `wFiltered_4DVolume.nii`, `w4DVolume.nii`, then `w4D.nii`. Motion files must contain one row per functional volume.

## Workflow

```matlab
addpath('matlab');
setup_study_paths;
cfg = study_config('/path/to/Data', '/path/to/derivatives');
cfg.spm_dir = '/path/to/spm12';
cfg.hcpex_source_atlas = '/path/to/HCPex_2mm.nii';

build_merged_hcpex_atlas(cfg);
cfg.hcpex_labels = fullfile(cfg.merged_atlas_dir, ...
    'HCPex_MergedROI_UrbanPerception_Refined17_lookup.csv');
cfg.hcpex_labels = update_atlas_voxel_counts(cfg);
cfg.subject_list = run_build_conditions(cfg);

run_activation_first_level(cfg);
run_activation_group_model(cfg);
run_activation_conjunction(cfg);

run_lss_first_level(cfg);
run_beta_series_connectivity(cfg);
Statistics = fit_connectivity_task_lme(cfg);
```

## Output layout

Activation products under `cfg.activation_root` include:

- `Atlas/`: nearest-neighbour functional-space atlas, ROI union mask and functional-grid ROI voxel counts.
- `Manifest/`: participant/run contrast images and model contrast definitions.
- `SecondLevel_TaskOneSample/Stim_gt_Fixation/<Task>/`: pooled one-sample task models.
- `SecondLevel_FullModel/Stim_gt_Fixation/`: repeated-measures model and contrasts.
- `FDR_Maps/`: statistic, p/q and significance-mask NIfTIs; positive task map stems are `Stim_gt_Fixation_OneSample_<Task>_gt_0`.
- `FDR_VoxelTables/` and `FDR_ROISummary/`: significant voxel lists and ROI summaries.
- `FourTask_ActivationConjunction/Stim_gt_Fixation/`: common activation, overlap labels and descriptive summaries.

LSS target images are stored at `cfg.lss_root/GroupN/RunN/UP_subject/LSS_StimBeta_RunN/beta_stim_trial_XXX.nii`. The event metadata associate `XXX` with its trial and image identity for the RSA bridge.

Connectivity products under `cfg.connectivity_root` include:

- `Step1_ROIBetaSeries/<Task>/*_betaSeries.mat`: trial-by-ROI matrix `X`, beta indices, ROI labels and source metadata.
- `Step2_BSCMatrices/<Task>/*_BSC.mat`: correlation `R`, Fisher-z `Z`, run metadata and invalid-ROI indicators.
- `BSC_Run_QC.csv`: beta counts and data-quality summaries.
- `BSC_EdgeLong_AllAvailable.csv`: edge values from available runs.
- `BSC_EdgeLong.csv`: complete-case participant/task/edge values.
- `BSC_all_tasks_stack_HCPex15.mat`: `Zall` indexed ROI-by-ROI-by-subject-by-task, with task names, subject keys and ROI metadata.
- `BSC_TaskMainEffect_LME.csv`: connection-level omnibus statistics, FDR-adjusted probabilities and model settings.
