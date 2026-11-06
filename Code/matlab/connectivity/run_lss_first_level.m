function run_lss_first_level(cfg)
% RUN_LSS_FIRST_LEVEL Estimate 99 target-trial beta images per retained run.
% Each model contains the target Stim, all other Stim trials, Rating, Fixation,
% Friston-24 motion terms, FD spike regressors and an intercept.
% The target-trial beta is collected in LSS_StimBeta_RunN for connectivity/RSA.
if nargin < 1, cfg = study_config(); end
required = {'spm_dir','subject_list','preprocessed_root','conditions_root','lss_root'};
for k = 1:numel(required)
    assert(isfield(cfg,required{k}) && ~isempty(cfg.(required{k})), 'Missing cfg.%s', required{k});
end
assert(isfile(cfg.subject_list),'Missing subject list: %s',cfg.subject_list);
assert(isfolder(cfg.preprocessed_root),'Missing preprocessed time-series root: %s',cfg.preprocessed_root);
assert(isfolder(cfg.conditions_root),'Build event files first: %s',cfg.conditions_root);
assert(isfolder(cfg.spm_dir),'SPM12 directory is unavailable: %s',cfg.spm_dir);
addpath(cfg.spm_dir);

spm('Defaults', 'fMRI');
spm_jobman('initcfg');

CFG0.SubjectListTxt = cfg.subject_list;
[CFG0.FirstLevelBase, CFG0.OutModelName] = fileparts(cfg.lss_root);
CFG0.FuncBase = cfg.preprocessed_root;
CFG0.MotionFolderName = 'RealignParameter';
CFG0.CondBase = cfg.conditions_root;

CFG0.FuncNiiCandidates = { ...
    'wFiltered_4DVolume.nii', ...
    'w4DVolume.nii', ...
    'w4D.nii' ...
    };

CFG0.TR      = 1;
CFG0.fmri_t  = 60;
CFG0.fmri_t0 = 30;

CFG0.HPF       = 128;
CFG0.UseDerivs = false;

CFG0.ModelTag = 'RT';
CFG0.ExpectedTrials = 99;
if isfield(cfg,'n_trials'), CFG0.ExpectedTrials = cfg.n_trials; end

CFG0.UseButtonReg = false;
CFG0.ButtonDur    = 0;         % 0 or 0.1

CFG0.UseFriston24 = true;
CFG0.UseFDScrub   = true;
CFG0.FD_Thr       = 0.5;
CFG0.FD_RadiusMM  = 50;        % convert rotation to displacement (approx.)

CFG0.SkipIfTrialEstimated = true;
CFG0.ReestimateTrial      = false;

CFG0.MakeStimContrast = true;  % create a simple Stim_i T-contrast
CFG0.CopyStimBeta     = true;  % collect beta image for Stim_i

if isfield(cfg,'lss_options')
    user = cfg.lss_options; names = fieldnames(user);
    for k = 1:numel(names)
        assert(isfield(CFG0,names{k}), 'Unknown lss_options option: %s',names{k});
        CFG0.(names{k}) = user.(names{k});
    end
end
PLAN = { ...
    'Group1', {'Run1','Run2','Run3','Run4'}; ...
    'Group2', {'Run1','Run2','Run3','Run4'}; ...
    'Group3', {'Run1','Run2','Run3','Run4'}; ...
    'Group4', {'Run1','Run2','Run3','Run4'}; ...
};
if isfield(cfg,'plan'), PLAN = cfg.plan; end

for gi = 1:size(PLAN,1)
    gName = PLAN{gi,1};
    runs  = PLAN{gi,2};

    for ri = 1:numel(runs)
        rName = runs{ri};

        CFG = CFG0;
        CFG.GroupName = gName;
        CFG.RunName   = rName;

        fprintf('\n====================================================\n');
        fprintf('[BATCH-LSS-LSS] %s | %s %s\n', CFG.OutModelName, CFG.GroupName, CFG.RunName);
        fprintf('====================================================\n');

        run_firstlevel_one_group_run_LSS(CFG);
    end
end

fprintf('\n[ALL DONE] LSS RT-clean first-level finished.\n');

end

function run_firstlevel_one_group_run_LSS(CFG)

    CFG.SubjectList = load_subjectlist_from_txt_exact(CFG.GroupName, CFG.RunName, CFG.SubjectListTxt);

    CFG.FirstLevelRoot = fullfile(CFG.FirstLevelBase, CFG.OutModelName, CFG.GroupName, CFG.RunName);
    if ~exist(CFG.FirstLevelRoot,'dir'); mkdir(CFG.FirstLevelRoot); end

    CFG.FuncFolderName = func_folder_by_run(CFG.RunName);
    CFG.FuncRoot   = fullfile(CFG.FuncBase, CFG.GroupName, CFG.RunName, CFG.FuncFolderName);
    CFG.MotionRoot = fullfile(CFG.FuncBase, CFG.GroupName, CFG.RunName, CFG.MotionFolderName);

    CFG.CondDir = fullfile(CFG.CondBase, sprintf('conds_glm_%s_%s_%s', ...
        CFG.GroupName, CFG.RunName, CFG.ModelTag));

    fail_log = {};
    subj_qc_rows = {};
    nOK = 0;
    nFail = 0;

    for i = 1:numel(CFG.SubjectList)
        sub = CFG.SubjectList{i};

        try
            subj_dir = fullfile(CFG.FirstLevelRoot, ['UP_' sub]);
            if ~exist(subj_dir,'dir'); mkdir(subj_dir); end

            func_dir = fullfile(CFG.FuncRoot, ['UP_' sub]);
            if ~exist(func_dir,'dir')
                error('Missing func dir: %s', func_dir);
            end

            nii4d = pick_first_existing(func_dir, CFG.FuncNiiCandidates);
            if isempty(nii4d)
                error('Missing 4D NIfTI in %s (checked: %s)', func_dir, strjoin(CFG.FuncNiiCandidates, ', '));
            end

            V = spm_vol(nii4d);
            nVol = numel(V);

            scans = cell(nVol,1);
            for v = 1:nVol
                scans{v} = sprintf('%s,%d', nii4d, v);
            end

            cond_file = fullfile(CFG.CondDir, sprintf('Condition_GLM_%s_%s_%s_%s.mat', ...
                CFG.ModelTag, CFG.GroupName, sub, CFG.RunName));

            if ~exist(cond_file,'file')
                error('Missing condition mat: %s', cond_file);
            end

            cond_spmready = make_cond_spm_ready(cond_file, subj_dir, nVol, CFG.TR);

            E = load(cond_spmready);
            names     = E.names;
            onsets    = E.onsets;
            durations = E.durations;

            if isstring(names), names = cellstr(names); end
            if ischar(names),   names = cellstr(names); end

            idxStim = find(strcmp(names,'Stim'),1);
            idxRating  = find(strcmp(names,'Rating'),1);
            idxFixation = find(strcmp(names,'Fixation'),1);

            if isempty(idxStim) || isempty(idxRating) || isempty(idxFixation)
                error('Condition mat must contain Stim + Rating + Fixation. Got: %s', strjoin(names, ', '));
            end

            stim_on  = onsets{idxStim}(:);
            stim_dur = durations{idxStim}(:);

            rating_on   = onsets{idxRating}(:);
            rating_dur  = durations{idxRating}(:);

            fixation_on  = onsets{idxFixation}(:);
            fixation_dur = durations{idxFixation}(:);

            nTrials = numel(stim_on);
            assert(nTrials==CFG.ExpectedTrials,'Expected %d trials in %s.',CFG.ExpectedTrials,cond_file);

            if isscalar(stim_dur), stim_dur = repmat(stim_dur, nTrials, 1); end
            if isscalar(rating_dur),  rating_dur  = repmat(rating_dur,  nTrials, 1); end
            if isscalar(fixation_dur), fixation_dur = repmat(fixation_dur, nTrials, 1); end

            valid_btn = isfinite(rating_on) & isfinite(rating_dur) & isfinite(fixation_on);
            button_on = rating_on(valid_btn) + rating_dur(valid_btn);
            button_dur = zeros(sum(valid_btn),1) + CFG.ButtonDur;

            reg_dir = fullfile(CFG.MotionRoot, ['UP_' sub]);
            if ~exist(reg_dir,'dir')
                error('Missing motion directory: %s', reg_dir);
            end

            task_id = sscanf(CFG.RunName, 'Run%d');
            rp_list = dir(fullfile(reg_dir, sprintf('rp_aUP_%s_BOLD_PA_TASK%d_*.txt', sub, task_id)));
            if isempty(rp_list)
                error('No rp file found for %s in %s', sub, reg_dir);
            end

            rp_path = '';
            for k = 1:numel(rp_list)
                cand = fullfile(reg_dir, rp_list(k).name);
                rp = load(cand);
                if size(rp,1) == nVol
                    rp_path = cand;
                    break;
                end
            end
            if isempty(rp_path)
                error('No rp file matches volumes for %s: volumes=%d', sub, nVol);
            end

            [nuis_path, FD, spike_idx, meanFD] = build_nuisance_file(rp_path, subj_dir, CFG);

            collect_dir = fullfile(subj_dir, sprintf('LSS_StimBeta_%s', CFG.RunName));
            if ~exist(collect_dir,'dir'); mkdir(collect_dir); end

            fprintf('[LSS-LSS] %s %s UP_%s | Trials=%d | meanFD=%.4f | nSpike=%d\n', ...
                CFG.GroupName, CFG.RunName, sub, nTrials, meanFD, numel(spike_idx));

            nDone = 0;

            for ti = 1:nTrials
                trial_dir = fullfile(subj_dir, sprintf('trial_%03d', ti));
                spmmat = fullfile(trial_dir, 'SPM.mat');

                if exist(spmmat,'file') && CFG.SkipIfTrialEstimated && ~CFG.ReestimateTrial
                    if CFG.CopyStimBeta
                        collect_stim_beta_if_exists(spmmat, collect_dir, ti);
                    end
                    if CFG.MakeStimContrast
                        make_contrast_stim_i(spmmat, ti);
                    end
                    continue;
                end

                if CFG.ReestimateTrial && exist(trial_dir,'dir')
                    try, rmdir(trial_dir,'s'); catch, end
                end
                if ~exist(trial_dir,'dir'); mkdir(trial_dir); end

                others = setdiff(1:nTrials, ti);

                clear jobs;
                cidx = 0;

                jobs{1}.spm.stats.fmri_spec.dir = {trial_dir};
                jobs{1}.spm.stats.fmri_spec.timing.units = 'secs';
                jobs{1}.spm.stats.fmri_spec.timing.RT = CFG.TR;
                jobs{1}.spm.stats.fmri_spec.timing.fmri_t  = CFG.fmri_t;
                jobs{1}.spm.stats.fmri_spec.timing.fmri_t0 = CFG.fmri_t0;

                jobs{1}.spm.stats.fmri_spec.sess.scans = scans;

                cidx = cidx + 1;
                jobs{1}.spm.stats.fmri_spec.sess.cond(cidx).name     = sprintf('Stim_%03d', ti);
                jobs{1}.spm.stats.fmri_spec.sess.cond(cidx).onset    = stim_on(ti);
                jobs{1}.spm.stats.fmri_spec.sess.cond(cidx).duration = stim_dur(ti);
                jobs{1}.spm.stats.fmri_spec.sess.cond(cidx).tmod     = 0;
                jobs{1}.spm.stats.fmri_spec.sess.cond(cidx).pmod     = struct('name', {}, 'param', {}, 'poly', {});
                jobs{1}.spm.stats.fmri_spec.sess.cond(cidx).orth     = 0;

                cidx = cidx + 1;
                jobs{1}.spm.stats.fmri_spec.sess.cond(cidx).name     = 'Stim_other';
                jobs{1}.spm.stats.fmri_spec.sess.cond(cidx).onset    = stim_on(others);
                jobs{1}.spm.stats.fmri_spec.sess.cond(cidx).duration = stim_dur(others);
                jobs{1}.spm.stats.fmri_spec.sess.cond(cidx).tmod     = 0;
                jobs{1}.spm.stats.fmri_spec.sess.cond(cidx).pmod     = struct('name', {}, 'param', {}, 'poly', {});
                jobs{1}.spm.stats.fmri_spec.sess.cond(cidx).orth     = 0;

                cidx = cidx + 1;
                jobs{1}.spm.stats.fmri_spec.sess.cond(cidx).name     = 'Rating';
                jobs{1}.spm.stats.fmri_spec.sess.cond(cidx).onset    = rating_on';
                jobs{1}.spm.stats.fmri_spec.sess.cond(cidx).duration = rating_dur';
                jobs{1}.spm.stats.fmri_spec.sess.cond(cidx).tmod     = 0;
                jobs{1}.spm.stats.fmri_spec.sess.cond(cidx).pmod     = struct('name', {}, 'param', {}, 'poly', {});
                jobs{1}.spm.stats.fmri_spec.sess.cond(cidx).orth     = 0;

                cidx = cidx + 1;
                jobs{1}.spm.stats.fmri_spec.sess.cond(cidx).name     = 'Fixation';
                jobs{1}.spm.stats.fmri_spec.sess.cond(cidx).onset    = fixation_on';
                jobs{1}.spm.stats.fmri_spec.sess.cond(cidx).duration = fixation_dur';
                jobs{1}.spm.stats.fmri_spec.sess.cond(cidx).tmod     = 0;
                jobs{1}.spm.stats.fmri_spec.sess.cond(cidx).pmod     = struct('name', {}, 'param', {}, 'poly', {});
                jobs{1}.spm.stats.fmri_spec.sess.cond(cidx).orth     = 0;

                if CFG.UseButtonReg
                    cidx = cidx + 1;
                    jobs{1}.spm.stats.fmri_spec.sess.cond(cidx).name     = 'Button';
                    jobs{1}.spm.stats.fmri_spec.sess.cond(cidx).onset    = button_on';
                    jobs{1}.spm.stats.fmri_spec.sess.cond(cidx).duration = button_dur';
                    jobs{1}.spm.stats.fmri_spec.sess.cond(cidx).tmod     = 0;
                    jobs{1}.spm.stats.fmri_spec.sess.cond(cidx).pmod     = struct('name', {}, 'param', {}, 'poly', {});
                    jobs{1}.spm.stats.fmri_spec.sess.cond(cidx).orth     = 0;
                end

                jobs{1}.spm.stats.fmri_spec.sess.multi = {''};
                jobs{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
                jobs{1}.spm.stats.fmri_spec.sess.multi_reg = {nuis_path};
                jobs{1}.spm.stats.fmri_spec.sess.hpf = CFG.HPF;

                jobs{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
                jobs{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
                if CFG.UseDerivs
                    jobs{1}.spm.stats.fmri_spec.bases.hrf.derivs = [1 0];
                end
                jobs{1}.spm.stats.fmri_spec.volt = 1;
                jobs{1}.spm.stats.fmri_spec.global = 'None';
                jobs{1}.spm.stats.fmri_spec.mthresh = 0.8;
                jobs{1}.spm.stats.fmri_spec.mask = {''};
                jobs{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

                jobs{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', ...
                    substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
                    substruct('.','spmmat'));
                jobs{2}.spm.stats.fmri_est.write_residuals = 0;
                jobs{2}.spm.stats.fmri_est.method.Classical = 1;

                spm_jobman('run', jobs);

                if CFG.CopyStimBeta
                    collect_stim_beta_if_exists(fullfile(trial_dir,'SPM.mat'), collect_dir, ti);
                end

                if CFG.MakeStimContrast
                    make_contrast_stim_i(fullfile(trial_dir,'SPM.mat'), ti);
                end

                nDone = nDone + 1;
                if mod(ti,10) == 0 || ti == nTrials
                    fprintf('  [LSS] %s %s UP_%s : trial %03d/%03d done\n', ...
                        CFG.GroupName, CFG.RunName, sub, ti, nTrials);
                end
            end

            fprintf('[OK] LSS RT-clean done: %s %s UP_%s | nTrials=%d | newly_estimated=%d\n', ...
                CFG.GroupName, CFG.RunName, sub, nTrials, nDone);

            subj_qc_rows(end+1,:) = { ...
                CFG.GroupName, CFG.RunName, sub, nVol, nTrials, ...
                meanFD, numel(spike_idx), CFG.FD_Thr, ...
                CFG.UseButtonReg, CFG.UseFriston24, CFG.UseFDScrub ...
                }; %#ok<AGROW>

            nOK = nOK + 1;

        catch ME
            fprintf('[FAIL] %s %s UP_%s | %s\n', CFG.GroupName, CFG.RunName, sub, ME.message);
            fail_log(end+1,:) = {CFG.GroupName, CFG.RunName, sub, ME.message}; %#ok<AGROW>
            nFail = nFail + 1;
        end
    end

    if ~isempty(fail_log)
        TF = cell2table(fail_log, 'VariableNames', {'Group','Run','Subject','Error'});
        out_fail = fullfile(CFG.FirstLevelRoot, sprintf('failed_subjects_LSS_%s_%s.csv', ...
            CFG.GroupName, CFG.RunName));
        writetable(TF, out_fail);
        fprintf('[WRITE] %s\n', out_fail);
    end

    if ~isempty(subj_qc_rows)
        nrow = size(subj_qc_rows, 1);
        
        Group         = strings(nrow,1);
        Run           = strings(nrow,1);
        Subject       = strings(nrow,1);
        NVol          = zeros(nrow,1);
        NTrial        = zeros(nrow,1);
        MeanFD        = zeros(nrow,1);
        N_FD_Spikes   = zeros(nrow,1);
        FD_Threshold  = zeros(nrow,1);
        UseButtonReg  = false(nrow,1);
        UseFriston24  = false(nrow,1);
        UseFDScrub    = false(nrow,1);
        
        for rr = 1:nrow
            Group(rr)        = string(subj_qc_rows{rr,1});
            Run(rr)          = string(subj_qc_rows{rr,2});
            Subject(rr)      = string(subj_qc_rows{rr,3});
            NVol(rr)         = double(subj_qc_rows{rr,4});
            NTrial(rr)       = double(subj_qc_rows{rr,5});
            MeanFD(rr)       = double(subj_qc_rows{rr,6});
            N_FD_Spikes(rr)  = double(subj_qc_rows{rr,7});
            FD_Threshold(rr) = double(subj_qc_rows{rr,8});
            UseButtonReg(rr) = logical(subj_qc_rows{rr,9});
            UseFriston24(rr) = logical(subj_qc_rows{rr,10});
            UseFDScrub(rr)   = logical(subj_qc_rows{rr,11});
        end
        
        TQ = table(Group, Run, Subject, NVol, NTrial, MeanFD, N_FD_Spikes, ...
            FD_Threshold, UseButtonReg, UseFriston24, UseFDScrub);
        
        out_qc = fullfile(char(CFG.FirstLevelRoot), ...
            sprintf('qc_LSS_%s_%s.csv', char(CFG.GroupName), char(CFG.RunName)));
        
        if iscell(out_qc)
            out_qc = out_qc{1};
        end
        out_qc = char(string(out_qc));
        
        fid = fopen(out_qc, 'w');
        if fid == -1
            error('Cannot open QC file for writing: %s', out_qc);
        end
        
        fprintf(fid, 'Group,Run,Subject,NVol,NTrial,MeanFD,N_FD_Spikes,FD_Threshold,UseButtonReg,UseFriston24,UseFDScrub\n');
        
        for rr = 1:size(subj_qc_rows,1)
            fprintf(fid, '%s,%s,%s,%d,%d,%.6f,%d,%.3f,%d,%d,%d\n', ...
                char(string(subj_qc_rows{rr,1})), ...
                char(string(subj_qc_rows{rr,2})), ...
                char(string(subj_qc_rows{rr,3})), ...
                double(subj_qc_rows{rr,4}), ...
                double(subj_qc_rows{rr,5}), ...
                double(subj_qc_rows{rr,6}), ...
                double(subj_qc_rows{rr,7}), ...
                double(subj_qc_rows{rr,8}), ...
                logical(subj_qc_rows{rr,9}), ...
                logical(subj_qc_rows{rr,10}), ...
                logical(subj_qc_rows{rr,11}));
        end
        
        fclose(fid);
        fprintf('[WRITE] %s\n', out_qc);
    end

    fprintf('[SUMMARY] LSS RT-clean %s %s | OK=%d | FAIL=%d\n', ...
        CFG.GroupName, CFG.RunName, nOK, nFail);
end

function collect_stim_beta_if_exists(spmmat, collect_dir, ti)
    try
        beta_path = find_beta_for_stim_i(spmmat, ti);
        if ~isempty(beta_path) && exist(beta_path,'file')
            dst = fullfile(collect_dir, sprintf('beta_stim_trial_%03d.nii', ti));
            copyfile(beta_path, dst);
        end
    catch
    end
end

function beta_path = find_beta_for_stim_i(spmmat, ti)
    beta_path = '';

    S = load(spmmat);
    SPM = S.SPM;

    xnames = SPM.xX.name;
    if isstring(xnames), xnames = cellstr(xnames); end
    if ischar(xnames),   xnames = cellstr(xnames); end

    xnorm = cellfun(@(s) normalize_xname(s), xnames, 'UniformOutput', false);

    target = sprintf('Stim_%03d', ti);
    pat = ['^Sn\(1\)\s*' regexptranslate('escape', target) '.*bf\(1\)$'];
    hit = ~cellfun('isempty', regexp(xnorm, pat, 'once'));
    idx = find(hit, 1);

    if isempty(idx)
        pat2 = ['^Sn\(1\)\s*' regexptranslate('escape', target)];
        hit2 = ~cellfun('isempty', regexp(xnorm, pat2, 'once'));
        idx = find(hit2, 1);
    end

    if isempty(idx)
        return;
    end

    beta_path = fullfile(fileparts(spmmat), sprintf('beta_%04d.nii', idx));
end

function make_contrast_stim_i(spmmat, ti)

    S = load(spmmat);
    SPM = S.SPM;

    xnames = SPM.xX.name;
    if isstring(xnames), xnames = cellstr(xnames); end
    if ischar(xnames),   xnames = cellstr(xnames); end

    xnorm = cellfun(@(s) normalize_xname(s), xnames, 'UniformOutput', false);

    target = sprintf('Stim_%03d', ti);
    pat = ['^Sn\(1\)\s*' regexptranslate('escape', target) '.*bf\(1\)$'];
    hit = ~cellfun('isempty', regexp(xnorm, pat, 'once'));
    idx = find(hit, 1);

    if isempty(idx)
        pat2 = ['^Sn\(1\)\s*' regexptranslate('escape', target)];
        hit2 = ~cellfun('isempty', regexp(xnorm, pat2, 'once'));
        idx = find(hit2, 1);
    end

    if isempty(idx)
        return;
    end

    ncol = size(SPM.xX.X,2);
    w = zeros(1, ncol);
    w(idx) = 1;

    SPM.xCon = [];
    SPM.xCon = spm_FcUtil('Set', sprintf('Stim_%03d', ti), 'T', 'c', w', SPM.xX.xKXs);
    spm_contrasts(SPM);
    save(spmmat, 'SPM');
end

function s = normalize_xname(s)
    if isstring(s), s = char(s); end
    if isempty(s), s = ''; return; end
    s = strrep(s, char(160), ' ');
    s = strtrim(s);
    s = regexprep(s, '\s+', ' ');
end

function [nuis_path, FD, spike_idx, meanFD] = build_nuisance_file(rp_path, out_dir, CFG)
% Expand the six realignment traces and their derivatives to 24 terms.
% FD spikes are one-hot regressors; rotations use a 50-mm displacement radius.

    rp = load(rp_path);
    if size(rp,2) ~= 6
        error('rp file must be nVol x 6: %s', rp_path);
    end

    nVol = size(rp,1);

    if CFG.UseFriston24
        d_rp = [zeros(1,6); diff(rp)];
        Xmot = [rp, d_rp, rp.^2, d_rp.^2];
    else
        Xmot = rp;
    end

    rot_mm = rp(:,4:6) * CFG.FD_RadiusMM;
    mot_mm = [rp(:,1:3), rot_mm];

    dmot = [zeros(1,6); diff(mot_mm)];
    FD = sum(abs(dmot), 2);
    meanFD = mean(FD);

    if CFG.UseFDScrub
        spike_idx = find(FD > CFG.FD_Thr);
        Xspk = zeros(nVol, numel(spike_idx));
        for k = 1:numel(spike_idx)
            Xspk(spike_idx(k), k) = 1;
        end
    else
        spike_idx = [];
        Xspk = zeros(nVol, 0);
    end

    Xnuis = [Xmot, Xspk];

    nuis_path = fullfile(out_dir, 'multi_reg_nuisance.txt');
    dlmwrite(nuis_path, Xnuis, 'delimiter', '\t', 'precision', '%.8f');

    fd_path = fullfile(out_dir, 'fd_trace.txt');
    dlmwrite(fd_path, FD, 'delimiter', '\t', 'precision', '%.8f');

    qc_path = fullfile(out_dir, 'nuisance_qc.mat');
    save(qc_path, 'rp', 'Xmot', 'FD', 'spike_idx', 'Xnuis');
end

function fp = pick_first_existing(dirpath, candidates)
    fp = '';
    for i = 1:numel(candidates)
        cand = fullfile(dirpath, candidates{i});
        if exist(cand,'file')
            fp = cand;
            return;
        end
    end
end

function folder = func_folder_by_run(runName)
    switch lower(runName)
        case 'run1'
            folder = 'FunImgARCFW';
        case 'run2'
            folder = 'S2_FunImgARCFW';
        case 'run3'
            folder = 'S3_FunImgARCFW';
        case 'run4'
            folder = 'S4_FunImgARCFW';
        otherwise
            error('Unknown RunName: %s', runName);
    end
end

function subjList = load_subjectlist_from_txt_exact(groupName, runName, txtPath)
    txtPath   = char(string(txtPath));
    groupName = char(string(groupName));
    runName   = char(string(runName));

    if ~exist(txtPath,'file')
        error('SubjectList.txt not found: %s', txtPath);
    end

    raw = fileread(txtPath);
    raw = strrep(raw, sprintf('\r\n'), sprintf('\n'));
    raw = strrep(raw, sprintf('\r'),   sprintf('\n'));

    lines = regexp(raw, '\n', 'split');
    if isstring(lines), lines = cellstr(lines); end
    lines = lines(:);

    function s = clean_line(s)
        if isstring(s), s = char(s); end
        s = strrep(s, char(65279), '');
        s = strrep(s, char(160), ' ');
        s = strrep(s, sprintf('\t'), ' ');
        s = strtrim(s);
    end

    clean = cell(size(lines));
    for i = 1:numel(lines)
        clean{i} = clean_line(lines{i});
    end

    gIdx = find(strcmpi(clean, groupName), 1, 'first');
    if isempty(gIdx)
        error('Cannot find block for %s in %s', groupName, txtPath);
    end

    endIdx = numel(clean) + 1;
    for i = gIdx+1:numel(clean)
        if ~isempty(regexp(clean{i}, '^Group\d+$', 'once'))
            endIdx = i;
            break;
        end
    end
    block = clean(gIdx+1:endIdx-1);

    runStart = [];
    for i = 1:numel(block)
        if ~isempty(regexp(block{i}, ['^' regexptranslate('escape', runName) '\s*='], 'once'))
            runStart = i;
            break;
        end
    end
    if isempty(runStart)
        error('Cannot find %s under %s in %s', runName, groupName, txtPath);
    end

    buf = block{runStart};
    j = runStart;
    while isempty(strfind(buf, '};')) && j < numel(block)
        j = j + 1;
        buf = [buf ' ' block{j}]; %#ok<AGROW>
    end

    tok = regexp(buf, ['^' regexptranslate('escape', runName) '\s*=\s*\{([^}]*)\}\s*;'], 'tokens', 'once');
    if isempty(tok) || isempty(tok{1})
        error('Found %s but cannot parse "{...}" list.', runName);
    end
    inside = tok{1};

    idTok = regexp(inside, '''(\d+)''', 'tokens');
    if isempty(idTok)
        error('Cannot parse quoted IDs under %s/%s.', groupName, runName);
    end

    subjList = cell(numel(idTok), 1);
    for k = 1:numel(idTok)
        sid = idTok{k}{1};
        n = str2double(sid);
        if isnan(n)
            error('Non-numeric subject ID: %s', sid);
        end
        subjList{k} = sprintf('%03d', n);
    end
end

function cond_out = make_cond_spm_ready(cond_in, out_dir, nVol, TR)

    S = load(cond_in);
    if ~isfield(S,'names') || ~isfield(S,'onsets') || ~isfield(S,'durations')
        error('Condition mat missing names/onsets/durations: %s', cond_in);
    end

    names = S.names;
    onsets = S.onsets;
    durations = S.durations;

    if isstring(names), names = cellstr(names); end
    if ischar(names),   names = cellstr(names); end
    if ~iscell(names)
        error('names must be cell.');
    end
    if ~iscell(onsets) || ~iscell(durations)
        error('onsets/durations must be cell.');
    end
    if numel(onsets) ~= numel(names) || numel(durations) ~= numel(names)
        error('names/onsets/durations length mismatch.');
    end

    t_scan = nVol * TR;
    t_cap  = t_scan - 1e-3;

    for c = 1:numel(names)
        o = double(onsets{c}(:));
        d = double(durations{c}(:));

        if numel(d) == 1 && numel(o) > 1
            d = repmat(d, size(o));
        end
        if numel(o) ~= numel(d)
            error('Cond %s: onset/duration mismatch.', names{c});
        end

        keep = isfinite(o) & isfinite(d);
        o = o(keep);
        d = d(keep);

        d(d < 0) = 0;

        keep2 = (o >= 0) & (o < t_scan);
        o = o(keep2);
        d = d(keep2);

        endt = o + d;
        over = endt > t_cap;
        if any(over)
            d(over) = max(0, t_cap - o(over));
        end

        onsets{c} = o(:)';
        durations{c} = d(:)';
    end

    cond_out = fullfile(out_dir, 'Condition_SPMREADY.mat');
    save(cond_out, 'names', 'onsets', 'durations');
end
