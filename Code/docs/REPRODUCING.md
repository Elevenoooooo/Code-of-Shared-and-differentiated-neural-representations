# Analysis workflow

## Configuration

From the Code directory, install the Python package with `python -m pip install -e .`. Configure MATLAB with:

```matlab
addpath('matlab'); setup_study_paths;
cfg = study_config('/study/Data', '/study/derivatives');
cfg.spm_dir = '/software/spm12';
cfg.dpabi_dir = '/software/DPABI';
```

Set cfg.dpabi_config_file to the DPABI Cfg MAT describing input sessions and acquisition-specific slice timing. Set the HCPex atlas, parcel crosswalk and feature-RDM locations in cfg. Outputs are organized under cfg.derivatives_root and cfg.output_root.

## 1. Preprocessing and ROI definition

Run `run_study_stage('preprocessing',cfg)` to apply initial-volume removal, slice timing, motion correction, tissue-based normalization, nuisance regression and filtering. The nuisance model includes white-matter, CSF and global signals, Friston-24 terms, high-motion indicators and a first-order trend. The mean signal is restored and the pass band is 0.01-0.10 Hz.

Run `cfg=run_study_stage('atlas',cfg)` to combine HCPex parcels into bilateral ROIs and write atlas metadata. Fifteen ROIs are retained for the main analyses. ROI labels are sampled using nearest-neighbour interpolation where image grids differ.

## 2. Events and conventional GLM

Run `cfg=run_study_stage('conditions',cfg)`. Each run retains 99 trials. Stimulus presentation is modelled for four seconds. The subsequent two-second interval is divided into the response period (Rating) and the remaining central-fixation period (Fixation). Trial order and behavioural response are parametric modulators on Stim; sequential pmod orthogonalization is disabled.

Run `run_study_stage('activation-first',cfg)` for first-level estimation with Friston-24 and motion-spike covariates. Stimulus > Fixation is the primary activation contrast; pmods and nuisance regressors have zero contrast weights.

## 3. Group activation and regional differentiation

Run `run_study_stage('activation-group',cfg)`. Positive task activation is evaluated with one-sample t tests. Repeated-measures task differences are represented with subject and group/task-cell effects. The models yield the omnibus task F statistic and six signed pairwise T contrasts. Voxelwise BH-FDR uses q<0.05.

Run `run_study_stage('conjunction',cfg)` to intersect the four task significance masks and calculate the minimum-T map. Activation-profile correlations use aligned unthresholded T values.

For each ROI and pairwise contrast, significance-weighted skewness is the FDR-significant voxel fraction multiplied by the adjusted skewness of the complete signed-T distribution. In VVC, PPC and vmPFC, the six SWS values are reconstructed into four zero-sum task coordinates by least squares. Euclidean distances summarize task separation.

```bash
python -m urban_perception.neural_summary sws --pairs /study/sws_pairs.csv --out outputs/sws
python -m urban_perception.neural_summary similarity --profiles /study/aligned_profiles.csv --out outputs/profile_similarity
```

The SWS table uses `TaskA,TaskB,ROI,SWS`. Profile tables use Safety,Vitality,Beauty,Imageability columns with aligned voxel or edge rows.

## 4. Beta-series connectivity

Run `run_study_stage('lss',cfg)` to estimate a target-stimulus beta for each of the 99 retained trials. Each LSS model includes the target, other stimuli, Rating, Fixation, motion covariates and an intercept.

Run `run_study_stage('connectivity',cfg)` to average beta values within the 15 ROIs, calculate Pearson correlations and Fisher-z transforms, and fit the task model for each of the 105 unique edges. The mixed-effects model is `Z ~ Task + Group + (1|Subject)`. BH-FDR is applied across the 105 omnibus task tests.

## 5. Representational similarity analysis

Prepare the ten environmental model RDMs:

```bash
python -m urban_perception.rsa --rdm-dir /study/feature_rdms --image-manifest config/image_id.csv --out-dir /study/derivatives/urban_rdms
```

Then run:

```matlab
run_study_stage('rsa-manifest', cfg);
run_study_stage('rsa-patterns', cfg);
run_study_stage('rsa', cfg);
```

Each participant/run/ROI yields a 99 x 99 neural RDM defined by 1 minus Pearson pattern correlation. Feature and rating RDMs use the same image order. Ten standardized feature-RDM vectors are fitted simultaneously in VVC, PPC and vmPFC. Rating RDMs are correlated with neural RDMs using Spearman correlation in all 15 retained ROIs. Participant coefficients/correlations enter group inference as described in [RSA methods](rsa.md).

## 6. Behaviour and gaze

Behavioural ratings and gaze maps are analysed with the Python modules described in [Behaviour and gaze methods](behavior_gaze_stimulus.md). The gaze pipeline applies calibration, valid-sample and duration criteria, constructs duration-weighted probability maps, and computes pooled Pearson and image-wise SSIM similarity. Confidence intervals and empirical tail probabilities use 5,000 image-level bootstrap samples; correction is separate across the six pairs for each metric.
