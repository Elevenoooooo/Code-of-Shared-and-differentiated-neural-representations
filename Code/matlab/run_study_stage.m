function cfg = run_study_stage(stage, cfg)
%RUN_STUDY_STAGE Run one stage of the study analysis workflow.
% Example: cfg = run_study_stage('conditions', cfg);
% The returned configuration carries a newly generated run-inclusion list.
if nargin < 2, cfg = study_config(); end
switch lower(char(stage))
    case 'preprocessing'
        run_dpabi_preprocessing(cfg);
    case 'atlas'
        build_merged_hcpex_atlas(cfg);
        cfg.hcpex_labels = fullfile(cfg.merged_atlas_dir, ...
            'HCPex_MergedROI_UrbanPerception_Refined17_lookup.csv');
        cfg.hcpex_labels = update_atlas_voxel_counts(cfg);
    case 'conditions'
        cfg.subject_list = run_build_conditions(cfg);
    case 'activation-first'
        run_activation_first_level(cfg);
    case 'activation-group'
        run_activation_group_model(cfg);
    case 'conjunction'
        run_activation_conjunction(cfg);
    case 'lss'
        run_lss_first_level(cfg);
    case 'connectivity'
        run_beta_series_connectivity(cfg);
        fit_connectivity_task_lme(cfg);
    case 'connectivity-lme'
        fit_connectivity_task_lme(cfg);
    case 'rsa-manifest'
        build_rsa_manifest(cfg);
    case 'rsa-patterns'
        extract_rsa_patterns(cfg);
    case 'rsa'
        run_rsa(cfg);
    otherwise
        error('Unknown stage: %s. See docs/REPRODUCING.md.', stage);
end
end
