function run_activation_first_level(cfg)
% RUN_ACTIVATION_FIRST_LEVEL Estimate stimulus-versus-fixation activation.
% The canonical-HRF model includes Stim, Rating, Fixation, two Stim pmods,
% Friston-24 realignment terms, FD spike regressors and an intercept.
% The primary T contrast assigns +1 to Stim, -1 to Fixation and zero elsewhere.
% Outputs: one SPM model, contrast images and motion summaries per subject/run.
if nargin < 1, cfg = study_config(); end
required = {'spm_dir','subject_list','preprocessed_root','conditions_root','first_level_root'};
for k = 1:numel(required)
    assert(isfield(cfg,required{k}) && ~isempty(cfg.(required{k})), 'Missing cfg.%s', required{k});
end
assert(isfile(cfg.subject_list),'Missing subject list: %s',cfg.subject_list);
assert(isfolder(cfg.preprocessed_root),'Missing preprocessed time-series root: %s',cfg.preprocessed_root);
assert(isfolder(cfg.conditions_root),'Build event files first: %s',cfg.conditions_root);
assert(isfolder(cfg.spm_dir),'SPM12 directory is unavailable: %s',cfg.spm_dir);
addpath(cfg.spm_dir);

spm('Defaults','fMRI');
spm_jobman('initcfg');

CFG0.SubjectListTxt = cfg.subject_list;
[CFG0.FirstLevelBase, CFG0.OutModelName] = fileparts(cfg.first_level_root);
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

CFG0.HPF        = 128;
CFG0.UseDerivs  = false;

CFG0.UseFriston24 = true;
CFG0.UseFDScrub   = true;
CFG0.FD_Thr       = 0.5;
CFG0.FD_RadiusMM  = 50;    % convert rotations (radians) to approximate displacement in mm

CFG0.ReestimateExistingModel = false;

CFG0.SkipIfAllContrastsExist = true;
CFG0.RecreateContrasts       = false;

MODEL_LIST = {'RT'};

CFG0.WantContrasts = {'Stim>Fixation'};

if isfield(cfg,'activation_first_level_options')
    user = cfg.activation_first_level_options; names = fieldnames(user);
    for k = 1:numel(names)
        assert(isfield(CFG0,names{k}), 'Unknown activation_first_level_options option: %s',names{k});
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

for mi = 1:numel(MODEL_LIST)
    modelTag = MODEL_LIST{mi};

    for gi = 1:size(PLAN,1)
        gName = PLAN{gi,1};
        runs  = PLAN{gi,2};

        for ri = 1:numel(runs)
            rName = runs{ri};

            CFG = CFG0;
            CFG.ModelTag  = modelTag;
            CFG.GroupName = gName;
            CFG.RunName   = rName;

            fprintf('\n====================================================\n');
            fprintf('[BATCH] Model=%s | %s %s\n', CFG.ModelTag, CFG.GroupName, CFG.RunName);
            fprintf('====================================================\n');

            run_firstlevel_one_group_run(CFG);
        end
    end
end

fprintf('\n[ALL DONE] First-level activation analysis finished.\n');

end

function run_firstlevel_one_group_run(CFG)

    CFG.SubjectList = load_subjectlist_from_txt_exact(CFG.GroupName, CFG.RunName, CFG.SubjectListTxt);

    if isfield(CFG, 'OutModelName') && ~isempty(CFG.OutModelName)
        outModelName = CFG.OutModelName;
    else
        outModelName = ['GLM_' CFG.ModelTag '_twoPmods'];
    end
    CFG.FirstLevelRoot = fullfile(CFG.FirstLevelBase, outModelName, CFG.GroupName, CFG.RunName);
    if ~exist(CFG.FirstLevelRoot,'dir'); mkdir(CFG.FirstLevelRoot); end

    CFG.FuncFolderName = func_folder_by_run(CFG.RunName);
    CFG.FuncRoot   = fullfile(CFG.FuncBase, CFG.GroupName, CFG.RunName, CFG.FuncFolderName);
    CFG.MotionRoot = fullfile(CFG.FuncBase, CFG.GroupName, CFG.RunName, CFG.MotionFolderName);

    CFG.CondDir = fullfile(CFG.CondBase, sprintf('conds_glm_%s_%s_%s', CFG.GroupName, CFG.RunName, CFG.ModelTag));

    fail_log = {};
    subj_qc_rows = {};
    nOK = 0; nFail = 0;

    for i = 1:numel(CFG.SubjectList)
        sub = CFG.SubjectList{i};

        try
            out_dir = fullfile(CFG.FirstLevelRoot, ['UP_' sub]);

            if isfield(CFG,'ReestimateExistingModel') && CFG.ReestimateExistingModel && exist(out_dir,'dir')
                try
                    rmdir(out_dir, 's');
                catch ME2
                    error('Cannot remove existing output dir: %s\n%s', out_dir, ME2.message);
                end
            end

            if ~exist(out_dir,'dir'); mkdir(out_dir); end

            spmmat = fullfile(out_dir, 'SPM.mat');

            if ~exist(spmmat,'file')
                func_dir = fullfile(CFG.FuncRoot, ['UP_' sub]);
                if ~exist(func_dir,'dir'); error('Missing func dir: %s', func_dir); end

                nii4d = pick_first_existing(fullfile(func_dir, ''), CFG.FuncNiiCandidates);
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
                sess_cond = load_condmat_to_sesscond(cond_file, nVol, CFG.TR);

                reg_dir = fullfile(CFG.MotionRoot, ['UP_' sub]);
                if ~exist(reg_dir,'dir'); error('Missing motion directory: %s', reg_dir); end

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

                clear jobs;
                jobs{1}.spm.stats.fmri_spec.dir = {out_dir};
                jobs{1}.spm.stats.fmri_spec.timing.units = 'secs';
                jobs{1}.spm.stats.fmri_spec.timing.RT = CFG.TR;
                jobs{1}.spm.stats.fmri_spec.timing.fmri_t  = CFG.fmri_t;
                jobs{1}.spm.stats.fmri_spec.timing.fmri_t0 = CFG.fmri_t0;

                [nuis_path, FD, spike_idx, meanFD] = build_nuisance_file(rp_path, out_dir, CFG);

                jobs{1}.spm.stats.fmri_spec.sess.scans  = scans;
                jobs{1}.spm.stats.fmri_spec.sess.cond  = sess_cond;
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
                    substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
                jobs{2}.spm.stats.fmri_est.write_residuals = 0;
                jobs{2}.spm.stats.fmri_est.method.Classical = 1;

                spm_jobman('run', jobs);
                fprintf('[EST] Model=%s %s %s UP_%s estimated | meanFD=%.4f | nSpike=%d | nuisance=%s\n', ...
                    CFG.ModelTag, CFG.GroupName, CFG.RunName, sub, meanFD, numel(spike_idx), nuis_path);
            else
                fprintf('[EST-SKIP] Model=%s %s %s UP_%s (SPM.mat exists)\n', CFG.ModelTag, CFG.GroupName, CFG.RunName, sub);
            end

            add_or_update_contrasts(spmmat, out_dir, CFG);

            qc_mat = fullfile(out_dir, 'nuisance_qc.mat');
            if exist(qc_mat, 'file')
                Q = load(qc_mat);
                if isfield(Q,'FD'), FD_for_qc = Q.FD; else, FD_for_qc = NaN; end
                if isfield(Q,'spike_idx'), spk_for_qc = Q.spike_idx; else, spk_for_qc = []; end
                meanFD_qc = mean_omitnan(FD_for_qc(:));
                nVol_qc = numel(FD_for_qc);
                nSpike_qc = numel(spk_for_qc);
            else
                meanFD_qc = NaN; nVol_qc = NaN; nSpike_qc = NaN;
            end

            subj_qc_rows(end+1,:) = {CFG.ModelTag, CFG.GroupName, CFG.RunName, sub, ...
                nVol_qc, meanFD_qc, nSpike_qc, CFG.FD_Thr, CFG.UseFriston24, CFG.UseFDScrub}; %#ok<AGROW>

            fprintf('[OK] Model=%s %s %s UP_%s\n', CFG.ModelTag, CFG.GroupName, CFG.RunName, sub);
            nOK = nOK + 1;

        catch ME
            fprintf('[FAIL] Model=%s %s %s UP_%s\n%s\n', CFG.ModelTag, CFG.GroupName, CFG.RunName, sub, getReport(ME,'extended','hyperlinks','off'));
            fail_log(end+1,:) = {CFG.ModelTag, CFG.GroupName, CFG.RunName, sub, ME.message}; %#ok<AGROW>
            nFail = nFail + 1;
        end
    end

    if ~isempty(fail_log)
        TF = cell2table(fail_log, 'VariableNames', {'Model','Group','Run','Subject','Error'});
        out_fail = fullfile(CFG.FirstLevelRoot, sprintf('failed_subjects_%s_%s_%s.csv', CFG.ModelTag, CFG.GroupName, CFG.RunName));
        writetable(TF, out_fail);
        fprintf('[WRITE] %s\n', out_fail);
    end

    if ~isempty(subj_qc_rows)
        TQ = cell2table(subj_qc_rows, 'VariableNames', ...
            {'Model','Group','Run','Subject','NVol','MeanFD','N_FD_Spikes','FD_Threshold','UseFriston24','UseFDScrub'});
        out_qc = fullfile(CFG.FirstLevelRoot, sprintf('qc_firstlevel_%s_%s_%s.csv', CFG.ModelTag, CFG.GroupName, CFG.RunName));
        writetable(TQ, out_qc);
        fprintf('[WRITE] %s\n', out_qc);
    end

    fprintf('[SUMMARY] Model=%s %s %s | OK=%d | FAIL=%d\n', CFG.ModelTag, CFG.GroupName, CFG.RunName, nOK, nFail);
end

function add_or_update_contrasts(spmmat, out_dir, CFG)

    S = load(spmmat);
    SPM = S.SPM;

    xnames = SPM.xX.name;
    if isstring(xnames), xnames = cellstr(xnames); end
    if ischar(xnames),   xnames = cellstr(xnames); end
    xnorm = cellfun(@(s) normalize_xname(s), xnames, 'UniformOutput', false);

    iStim  = find_condition_main_effect(xnorm, 'Stim', {'order','respbin'});
    iFixation  = find_condition_main_effect(xnorm, 'Fixation', {});
    iOrder = find_pmod_effect(xnorm, 'Stim', 'Order');
    iResp  = find_pmod_effect(xnorm, 'Stim', 'RespBin');

    dump = fullfile(out_dir, 'design_xX_name_dump.txt');

    if isempty(iStim) || isempty(iFixation) || isempty(iOrder) || isempty(iResp)
        fid = fopen(dump, 'w');
        if fid > 0
            for k = 1:min(numel(xnorm), 1000)
                fprintf(fid, '%s\n', xnorm{k});
            end
            fclose(fid);
        end
        error(['Cannot locate regressors. ', ...
               'Stim=%d, Fixation=%d, Order=%d, RespBin=%d. ', ...
               'Dumped xX.name to: %s'], ...
               numel(iStim), numel(iFixation), numel(iOrder), numel(iResp), dump);
    end

    iStim  = iStim(1);
    iFixation  = iFixation(1);
    iOrder = iOrder(1);
    iResp  = iResp(1);

    ncol = size(SPM.xX.X,2);

    C = struct([]);
    C(1).name = 'Stim>Fixation'; C(1).w = make_weight(ncol, [iStim, iFixation], [1, -1]);

    existNames = {};
    if isfield(SPM,'xCon') && ~isempty(SPM.xCon)
        existNames = {SPM.xCon.name};
    end

    if CFG.SkipIfAllContrastsExist && ~CFG.RecreateContrasts
        allExist = all(ismember(CFG.WantContrasts, existNames));
        if allExist
            fprintf('[CON-SKIP] all desired contrasts already exist in %s\n', out_dir);
            return;
        end
    end

    clear conjob;
    conjob{1}.spm.stats.con.spmmat = {spmmat};
    cc = 0;
    for i = 1:numel(C)
        if ~ismember(C(i).name, CFG.WantContrasts)
            continue;
        end
        cc = cc + 1;
        conjob{1}.spm.stats.con.consess{cc}.tcon.name    = C(i).name;
        conjob{1}.spm.stats.con.consess{cc}.tcon.weights = C(i).w;
        conjob{1}.spm.stats.con.consess{cc}.tcon.sessrep = 'none';
    end
    conjob{1}.spm.stats.con.delete = CFG.RecreateContrasts;

    spm_jobman('run', conjob);

    fprintf('[CON] %s ensured:\n', out_dir);
    fprintf('      Stim=%d, Fixation=%d, Order=%d, RespBin=%d\n', iStim, iFixation, iOrder, iResp);
end

function idx = find_condition_main_effect(xnorm, condName, excludeWords)
    idx = [];
    condName2 = regexprep(lower(condName),'_','[_\s]*');

    cands = [];
    for k = 1:numel(xnorm)
        s = lower(xnorm{k});
        if ~isempty(regexp(s, ['^sn\(\d+\)\s+' condName2 '.*bf\(1\)$'], 'once'))
            cands(end+1) = k; %#ok<AGROW>
        end
    end
    if isempty(cands), return; end

    keep = true(size(cands));
    for i = 1:numel(cands)
        s = lower(xnorm{cands(i)});
        for w = 1:numel(excludeWords)
            if contains(s, lower(excludeWords{w}))
                keep(i) = false;
            end
        end
        if contains(s, 'x') && (contains(s,'order') || contains(s,'respbin') || contains(s,'^'))
            keep(i) = false;
        end
    end
    idx = cands(keep);
end

function idx = find_pmod_effect(xnorm, condName, modName)
    idx = [];
    cands = [];
    for k = 1:numel(xnorm)
        s = lower(xnorm{k});
        if contains(s, lower(condName)) && contains(s, lower(modName)) && contains(s, 'bf(1)')
            cands(end+1) = k; %#ok<AGROW>
        end
    end
    idx = cands;
end

function w = make_weight(ncol, idx, val)
    w = zeros(1, ncol);
    if numel(idx) ~= numel(val)
        if isscalar(val)
            val = repmat(val, size(idx));
        else
            error('Weight length mismatch.');
        end
    end
    for i = 1:numel(idx)
        w(idx(i)) = val(i);
    end
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

    if isfield(CFG,'UseFriston24') && CFG.UseFriston24
        d_rp = [zeros(1,6); diff(rp)];
        Xmot = [rp, d_rp, rp.^2, d_rp.^2];
    else
        Xmot = rp;
    end

    if ~isfield(CFG,'FD_RadiusMM') || isempty(CFG.FD_RadiusMM)
        CFG.FD_RadiusMM = 50;
    end
    rot_mm = rp(:,4:6) * CFG.FD_RadiusMM;
    mot_mm = [rp(:,1:3), rot_mm];

    dmot = [zeros(1,6); diff(mot_mm)];
    FD = sum(abs(dmot), 2);
    meanFD = mean(FD);

    if isfield(CFG,'UseFDScrub') && CFG.UseFDScrub
        if ~isfield(CFG,'FD_Thr') || isempty(CFG.FD_Thr)
            CFG.FD_Thr = 0.5;
        end
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
    save(qc_path, 'rp', 'Xmot', 'FD', 'spike_idx', 'Xnuis', 'rp_path', 'nuis_path');
end

function m = mean_omitnan(x)
    x = double(x(:));
    x = x(isfinite(x));
    if isempty(x)
        m = NaN;
    else
        m = mean(x);
    end
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
        case 'run1', folder = 'FunImgARCFW';
        case 'run2', folder = 'S2_FunImgARCFW';
        case 'run3', folder = 'S3_FunImgARCFW';
        case 'run4', folder = 'S4_FunImgARCFW';
        otherwise, error('Unknown RunName: %s', runName);
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
            endIdx = i; break;
        end
    end
    block = clean(gIdx+1:endIdx-1);

    runStart = [];
    for i = 1:numel(block)
        if ~isempty(regexp(block{i}, ['^' regexptranslate('escape', runName) '\s*='], 'once'))
            runStart = i; break;
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
        if isnan(n), error('Non-numeric subject ID: %s', sid); end
        subjList{k} = sprintf('%03d', n);
    end
end

function sess_cond = load_condmat_to_sesscond(cond_in, nVol, TR)

    S = load(cond_in);
    if ~isfield(S,'names') || ~isfield(S,'onsets') || ~isfield(S,'durations')
        error('Condition mat missing names/onsets/durations: %s', cond_in);
    end

    names = S.names;
    onsets = S.onsets;
    durations = S.durations;

    if isstring(names), names = cellstr(names); end
    if ischar(names),   names = cellstr(names); end
    if ~iscell(names),  error('names must be cellstr: %s', cond_in); end

    if ~iscell(onsets) || ~iscell(durations)
        error('onsets/durations must be cell arrays: %s', cond_in);
    end

    nCond = numel(names);
    if numel(onsets)~=nCond || numel(durations)~=nCond
        error('names/onsets/durations length mismatch: %s', cond_in);
    end

    if isfield(S,'orth')
        orth = S.orth;
        if iscell(orth)
            orth = double(cell2mat(orth(:)))';
        else
            orth = double(orth(:))';
        end
        if isempty(orth), orth = ones(1,nCond); end
        if numel(orth)~=nCond, orth = repmat(orth(1),1,nCond); end
    else
        orth = ones(1,nCond);
    end

    t_scan = nVol * TR;
    t_cap  = t_scan - 1e-3;

    empty_pmod = struct('name', {}, 'param', {}, 'poly', {});
    sess_cond  = repmat(struct('name','','onset',[],'duration',[],'tmod',0,'pmod',empty_pmod,'orth',1), 1, nCond);

    hasPmod = isfield(S,'pmod');

    for c = 1:nCond
        o0 = double(onsets{c}(:));
        d0 = double(durations{c}(:));
        if numel(d0)==1 && numel(o0)>1, d0 = repmat(d0,size(o0)); end
        if numel(o0)~=numel(d0)
            error('Onset/duration mismatch in %s (%s)', cond_in, names{c});
        end

        keep = isfinite(o0) & isfinite(d0);
        o = o0(keep); d = d0(keep);
        d(d<0) = 0;

        keep2 = (o>=0) & (o<t_scan);
        o = o(keep2); d = d(keep2);

        endt = o + d;
        over = endt > t_cap;
        if any(over), d(over) = max(0, t_cap - o(over)); end

        keep3 = (d > 0) | (o < t_cap);
        o = o(keep3); d = d(keep3);

        idx1 = find(keep);
        idx2 = idx1(keep2);
        idx3 = idx2(keep3);

        sess_cond(c).name     = char(string(names{c}));
        sess_cond(c).onset    = o(:)';
        sess_cond(c).duration = d(:)';
        sess_cond(c).tmod     = 0;
        sess_cond(c).orth     = double(orth(c) ~= 0);
        sess_cond(c).pmod     = empty_pmod;

        if hasPmod && c <= numel(S.pmod)
            pm = S.pmod(c);

            if isfield(pm,'name') && ~isempty(pm.name)
                if iscell(pm.name)
                    pm_names = pm.name;
                else
                    pm_names = {char(string(pm.name))};
                end

                if iscell(pm.param)
                    pm_params = pm.param;
                else
                    pm_params = {double(pm.param(:)')};
                end

                if iscell(pm.poly)
                    pm_polys = cellfun(@(x) double(x(1)), pm.poly);
                else
                    pm_polys = double(pm.poly(1));
                end
                if isscalar(pm_polys)
                    pm_polys = repmat(pm_polys, 1, numel(pm_names));
                end

                pmods = repmat(struct('name','','param',[],'poly',1), 1, numel(pm_names));

                for ii = 1:numel(pm_names)
                    pmods(ii).name = char(string(pm_names{ii}));
                    pk = double(pm_params{ii}(:));
                    if numel(pk) == numel(o0)
                        pk = pk(idx3);
                    end
                    pmods(ii).param = pk(:);
                    pmods(ii).poly  = double(pm_polys(ii));
                end

                sess_cond(c).pmod = pmods;
            end
        end
    end

    fprintf('[COND->SESSCOND] Built %d conditions from %s\n', nCond, cond_in);
end
