# Urban visual perception

Analysis code for **Shared and differentiated neural representations of urban scenes in the human brain**.

The study examines safety, vitality, beauty and imageability using behavioural ratings, eye tracking and fMRI. The workflow covers image-level response profiles, gaze-map similarity, task activation and conjunction, beta-series functional connectivity, significance-weighted skewness, and representational similarity analysis.

## Repository structure

```text
Code/
  config/                  study parameters, participant and image manifests, ROI definitions
  matlab/
    preprocessing/         DPABI preprocessing
    atlas/                 bilateral HCPex ROI construction
    activation/            events, first-level GLM, group statistics and conjunction
    connectivity/          LSS estimation, beta-series correlations and mixed-effects models
    rsa/                   trial-pattern extraction and run-level RSA
  python/urban_perception/ behavioural, gaze, image-feature and neural-summary analyses
  docs/                    methods, data formats and workflow documentation
  CITATION.cff             citation metadata
  LICENSE                  MIT license for the custom analysis code
```

## Software requirements

- Python 3.10 or later with NumPy, pandas, SciPy, scikit-image, Pillow and Matplotlib.
- MATLAB R2022b with Statistics and Machine Learning Toolbox and applicable image-processing functions.
- SPM12 for fMRI model estimation and NIfTI operations.
- DPABI for fMRI preprocessing.
- HCPex atlas files for ROI construction.

Install the Python package from this directory:

```bash
python -m pip install -e .
```

For Python NIfTI summaries, install `python -m pip install -e ".[neuro]"`. For the optional image-feature utilities, install `python -m pip install -e ".[features]"`. Dependency ranges are specified in `pyproject.toml` and `requirements.txt`. A Conda environment is provided in `environment.yml`.

## Study configuration

The companion data directory contains `Images`, `Behavior`, `ETHeatmap` and `fMRI`. Configure additional analysis inputs and software installations with `config/study_config.m`.

```matlab
addpath('matlab');
setup_study_paths;
cfg = study_config('/study/Data', '/study/derivatives');
cfg.spm_dir = '/software/spm12';
cfg.dpabi_dir = '/software/DPABI';
```

`config/participants.csv` identifies participants by Group and Subject. `config/subjects.txt` supplies the corresponding group/run lists for MATLAB. `config/image_id.csv` contains stimulus identifiers and experimental indices. `config/roi_definitions.csv` describes the 17 bilateral ROI labels and the 15 ROIs retained for analysis.

## Analysis workflow

| Stage | Implementation | Main output |
|---|---|---|
| Behavioural evaluation | `urban_perception.behavior` | Image-level positive-response rates and six cross-task correlations |
| Visual attention | `urban_perception.gaze` | Duration-weighted probability maps, Pearson/SSIM similarity and bootstrap inference |
| Preprocessing and ROIs | `run_dpabi_preprocessing`, `build_merged_hcpex_atlas` | Preprocessed time series and bilateral ROI masks |
| Task activation | `run_build_conditions`, `run_activation_first_level`, `run_activation_group_model` | Stimulus > Fixation contrasts, task T maps, omnibus F maps and pairwise differences |
| Common activation | `run_activation_conjunction` | Four-task conjunction and minimum-T map |
| Functional connectivity | `run_lss_first_level`, `run_beta_series_connectivity`, `fit_connectivity_task_lme` | Fisher-z connectivity for 105 edges and task-effect inference |
| Regional differentiation | `urban_perception.neural_summary` | SWS, zero-sum task coordinates and distances |
| Representational structure | `build_rsa_manifest`, `extract_rsa_patterns`, `run_rsa` | 99 x 99 neural RDMs, environmental-feature RSA-GLM and rating RSA |

### Behaviour

```bash
python -m urban_perception.behavior binary --behavior-dir ../Data/Behavior --participants config/participants.csv --output-dir outputs/behavior
python -m urban_perception.behavior online --ratings /study/online_image_means.csv --image-column File_Name --output-dir outputs/online
```

### Gaze

```bash
python -m urban_perception.gaze from-fixations --trials /study/trials_qc.csv --fixations /study/fixations.csv --participants config/participants.csv --output-dir outputs/gaze_maps
python -m urban_perception.gaze similarity --index /study/group_map_index.csv --n-bootstrap 5000 --output-dir outputs/gaze_similarity
```

### MATLAB stages

```matlab
run_study_stage('preprocessing', cfg);
cfg = run_study_stage('atlas', cfg);
cfg = run_study_stage('conditions', cfg);
run_study_stage('activation-first', cfg);
run_study_stage('activation-group', cfg);
run_study_stage('conjunction', cfg);
run_study_stage('lss', cfg);
run_study_stage('connectivity', cfg);
run_study_stage('rsa-manifest', cfg);
run_study_stage('rsa-patterns', cfg);
run_study_stage('rsa', cfg);
```

Prepare environmental feature RDMs before the RSA stage:

```bash
python -m urban_perception.rsa --rdm-dir /study/feature_rdms --image-manifest config/image_id.csv --out-dir /study/derivatives/urban_rdms
```

The 99 retained trials in each run define a 99 x 99 neural RDM. Run-specific feature and rating RDMs use the same stimulus order.

## Documentation

- [Workflow](docs/REPRODUCING.md)
- [Data formats and configuration](docs/DATA_REQUIREMENTS.md)
- [Behaviour, gaze and stimulus methods](docs/behavior_gaze_stimulus.md)
- [Activation and connectivity methods](docs/activation_connectivity.md)
- [Representational similarity analysis](docs/rsa.md)
- [Figure-to-code index](docs/FIGURE_CODE_MAP.md)

## Citation and license

Citation metadata are supplied in `CITATION.cff`. Custom analysis code is released under the MIT license. Software packages, atlases and other external resources retain their respective licenses; see [Third-party resources](docs/THIRD_PARTY.md).
