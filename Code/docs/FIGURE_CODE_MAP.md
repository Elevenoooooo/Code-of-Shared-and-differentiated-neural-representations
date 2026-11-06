# Figure-to-code index

| Manuscript analysis | Implementation | Output |
|---|---|---|
| Fig. 1: stimulus characteristics | `stimulus.py`, `features.py` | Feature-distribution and image-descriptor summaries |
| Fig. 2b: cross-dimensional ratings | `behavior.py` | Image-level rates and pairwise Pearson correlations |
| Fig. 2c: visual attention | `gaze.py` | Probability maps, Pearson/SSIM, bootstrap and BH results |
| Fig. 3a: task and shared activation | `run_activation_first_level.m`, `run_activation_group_model.m`, `run_activation_conjunction.m` | Task T maps, FDR masks and minimum-T conjunction |
| Fig. 3b: shared connectivity | `run_lss_first_level.m`, `run_beta_series_connectivity.m` | 15-ROI beta-series connectivity matrices |
| Fig. 3c: neural-profile similarity | `neural_summary.py` | Six activation- or connectivity-profile correlations |
| Fig. 4: activation differentiation | `run_activation_group_model.m`, `neural_summary.py` | Omnibus F, pairwise signed T, SWS and task coordinates |
| Fig. 4: connectivity differentiation | `fit_connectivity_task_lme.m` | Edge-level task F tests and BH-adjusted results |
| Fig. 5: environmental-feature RSA | `rsa.py`, `build_rsa_manifest.m`, `extract_rsa_patterns.m`, `run_rsa.m` | 99 x 99 run RDMs and simultaneous ten-feature model coefficients |
| Fig. 5: rating RSA | `run_rsa.m`, `rsa.py` | Task-specific neural-rating Spearman associations |

The analyses produce statistical maps, model tables and summary coordinates for figure composition.
