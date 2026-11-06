function run_dpabi_preprocessing(cfg)
%RUN_DPABI_PREPROCESSING Apply the study preprocessing settings with DPABI.
% The Cfg MAT supplies acquisition-specific timing, tissue masks and sessions.
assert(isfolder(cfg.dpabi_dir), 'Set cfg.dpabi_dir or DPABI_DIR.');
assert(isfolder(cfg.spm_dir), 'Set cfg.spm_dir or SPM_DIR.');
assert(isfile(cfg.dpabi_config_file), 'Provide cfg.dpabi_config_file.');
addpath(cfg.spm_dir);
addpath(genpath(cfg.dpabi_dir));
input = load(cfg.dpabi_config_file, 'Cfg');
Cfg = input.Cfg;
Cfg.TR = cfg.TR;
Cfg.IsRemoveFirstTimePoints = 1;
Cfg.RemoveFirstTimePoints = 6;
Cfg.IsSliceTiming = 1;
Cfg.SliceTiming.SliceNumber = cfg.n_slices;
if ~isempty(cfg.slice_order), Cfg.SliceTiming.SliceOrder = cfg.slice_order; end
assert(numel(Cfg.SliceTiming.SliceOrder) == cfg.n_slices, ...
    'Slice timing must specify all acquired slices.');
Cfg.SliceTiming.ReferenceSlice = cfg.reference_slice;
Cfg.IsRealign = 1;
Cfg.IsNeedT1CoregisterToFun = 1;
Cfg.IsSegment = 2;
Cfg.IsDARTEL = 1;
Cfg.IsNormalize = 3;
Cfg.Normalize.VoxSize = [2.4 2.4 2.4];
Cfg.IsCovremove = 1;
Cfg.Covremove.Timing = 'AfterNormalize';
Cfg.Covremove.PolynomialTrend = 1;
Cfg.Covremove.HeadMotion = 4; % Friston-24 expansion.
Cfg.Covremove.IsHeadMotionScrubbingRegressors = 1;
Cfg.Covremove.HeadMotionScrubbingRegressors.FDThreshold = cfg.fd_threshold;
Cfg.Covremove.WM.IsRemove = 1;
Cfg.Covremove.CSF.IsRemove = 1;
Cfg.Covremove.WholeBrain.IsRemove = 1;
Cfg.Covremove.WholeBrain.IsBothWithWithoutGSR = 0;
Cfg.Covremove.IsAddMeanBack = 1;
Cfg.IsFilter = 1;
Cfg.Filter.Timing = 'AfterNormalize';
Cfg.Filter.AHighPass_LowCutoff = 0.01;
Cfg.Filter.ALowPass_HighCutoff = 0.10;
Cfg.Filter.AAddMeanBack = 'Yes';
DPARSFA_run(Cfg);
end
