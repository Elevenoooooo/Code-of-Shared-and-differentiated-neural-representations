function generatedSubjectList = run_build_conditions(cfg)
% RUN_BUILD_CONDITIONS Create stimulus, rating and fixation events for SPM.
% Each run retains 99 trials after removal of the first six-second trial.
% Rating lasts until the response; Fixation fills the remaining rating interval.
% Stim carries standardized trial-order and signed-response pmods (orth=0).
% Returns the generated valid-response subject-list path. Trial metadata link
% each beta index to the observed image and original behavioural trial.
if nargin < 1, cfg = study_config(); end
required = {'subject_list','data_root','conditions_root','image_id_csv'};
for k = 1:numel(required)
    assert(isfield(cfg,required{k}) && ~isempty(cfg.(required{k})), 'Missing cfg.%s',required{k});
end
if ~isfolder(cfg.conditions_root), mkdir(cfg.conditions_root); end

CFG.SubjectListTxt = cfg.subject_list;
CFG.BehaviorRoot = fullfile(cfg.data_root,'Behavior');
CFG.OutBase = cfg.conditions_root;

CFG.RetainedTrials = 99;
if isfield(cfg,'n_trials'), CFG.RetainedTrials = cfg.n_trials; end
CFG.nTrials = CFG.RetainedTrials+1;
CFG.trialLen = 6;      % 4s stim + 2s rating
CFG.stimDur  = 4;
CFG.rateDur  = 2;

CFG.TR = 1;
CFG.DropFirstVols = 6; % first 6s removed in NIfTI
CFG.DropTrials = round((CFG.DropFirstVols * CFG.TR) / CFG.trialLen); % should be 1

CFG.ModelTag = 'RT';
CFG.IncludeRating = true;

CFG.ImageID_CSV = cfg.image_id_csv;
CFG.ImageID_FileColCandidates = {'File_Name','file_name','filename','File','Path','Path_raw'};
CFG.ImageID_IDColCandidates   = {'Img_ID','img_id','ID','ImageID','image_id'};

CFG.AddStimOrderPmod = true;
CFG.OrderZScore      = true;
CFG.OrderPoly        = 1;

CFG.AddStimRespPmod  = true;
CFG.RespPoly         = 1;
CFG.RespMissingValue = 0;  % missing/invalid response -> 0
CFG.StimOrth         = 0;  % IMPORTANT: disable sequential orth on Stim pmods

CFG.WriteStimMetaCSV = true;

CFG.RT_ColCandidates = { ...
    'key_resp_sc.rt_raw', ...
    'key_resp_sc_rt_raw', ...
    'key_resp_sc_rt', ...
    'rt_raw', ...
    'rt' ...
};

CFG.KeyRespColCandidates = { ...
    'key_resp_sc.keys_raw', ...
    'key_resp_sc_keys_raw', ...
    'key_resp_sc.keys', ...
    'key_resp_sc_keys', ...
    'keys_raw', ...
    'keys' ...
};

CFG.RT_Min = 0;
CFG.RT_Max = 2;

CFG.UpdatedSubjectListTxt = fullfile(cfg.conditions_root,'subject_list_generated_valid_keys.txt');
CFG.ExclusionLogCSV       = fullfile(CFG.OutBase, 'conditions_excluded_runs.csv');
CFG.QCSummaryCSV          = fullfile(CFG.OutBase, 'conditions_summary.csv');

if isfield(cfg,'condition_options')
    fields = fieldnames(cfg.condition_options);
    for k = 1:numel(fields)
        assert(isfield(CFG,fields{k}),'Unknown condition_options key: %s',fields{k});
        CFG.(fields{k}) = cfg.condition_options.(fields{k});
    end
    if ~isfield(cfg.condition_options,'DropTrials')
        CFG.DropTrials = round(CFG.DropFirstVols*CFG.TR/CFG.trialLen);
    end
end
assert(~strcmpi(char(CFG.UpdatedSubjectListTxt),char(CFG.SubjectListTxt)), ...
    'Generated subject list must not overwrite the bundled input list.');
assert(isfolder(CFG.BehaviorRoot),'Missing Behavior directory: %s',CFG.BehaviorRoot);
assert(mod(CFG.DropFirstVols*CFG.TR,CFG.trialLen)==0, ...
    'Discarded scan duration must equal whole trials.');
generatedSubjectList = CFG.UpdatedSubjectListTxt;

PLAN = { ...
    'Group1', {'Run1','Run2','Run3','Run4'}; ...
    'Group2', {'Run1','Run2','Run3','Run4'}; ...
    'Group3', {'Run1','Run2','Run3','Run4'}; ...
    'Group4', {'Run1','Run2','Run3','Run4'}; ...
};
if isfield(cfg,'plan'), PLAN = cfg.plan; end
TaskTagMap = tasktag_map();

SubjectMap = parse_subjectlist_txt(CFG.SubjectListTxt);
IMG = load_imageid_maps(CFG.ImageID_CSV, CFG.ImageID_FileColCandidates, CFG.ImageID_IDColCandidates);
fprintf('[IMG] Loaded imageID map: N=%d | fileCol=%s | idCol=%s\n', IMG.N, IMG.fileCol, IMG.idCol);

fprintf('\n====================================================\n');
fprintf('[STEP 1] Filtering runs with NO key_resp_sc.keys_raw in usable trials\n');
fprintf('====================================================\n');

FilteredMap = SubjectMap;
excl_rows = {};

for gi = 1:size(PLAN,1)
    gName = PLAN{gi,1};
    runs  = PLAN{gi,2};

    for ri = 1:numel(runs)
        rName = runs{ri};
        taskTag = TaskTagMap.(gName).(rName);

        if ~isfield(SubjectMap, gName) || ~isfield(SubjectMap.(gName), rName)
            warning('Missing %s/%s in subject list.', gName, rName);
            FilteredMap.(gName).(rName) = {};
            continue;
        end

        subjList = SubjectMap.(gName).(rName);
        keepSubs = {};
        behDir = fullfile(CFG.BehaviorRoot, gName, rName);

        for si = 1:numel(subjList)
            sub = subjList{si};
            csv_file = fullfile(behDir, sprintf('%s_%s_%s_%s_trials.csv', ...
                gName, sub, lower(rName), taskTag));

            if ~exist(csv_file, 'file')
                excl_rows(end+1,:) = {gName,rName,sub,taskTag,'CSV_MISSING'}; %#ok<SAGROW>
                continue;
            end

            try
                T = readtable(csv_file);
            catch ME
                excl_rows(end+1,:) = {gName,rName,sub,taskTag,['READ_FAIL: ' ME.message]}; %#ok<SAGROW>
                continue;
            end

            if ~any(strcmp(T.Properties.VariableNames,'thisTrialN_raw'))
                excl_rows(end+1,:) = {gName,rName,sub,taskTag,'MISSING_thisTrialN_raw'}; %#ok<SAGROW>
                continue;
            end

            [hasKeyCol, keyColName] = find_first_existing_var(T, CFG.KeyRespColCandidates);
            if ~hasKeyCol
                excl_rows(end+1,:) = {gName,rName,sub,taskTag,'MISSING_key_resp_sc.keys_raw'}; %#ok<SAGROW>
                continue;
            end

            use_mask = (T.thisTrialN_raw >= 0) & (T.thisTrialN_raw <= CFG.nTrials-1);
            Tuse = T(use_mask,:);
            if height(Tuse) ~= CFG.nTrials
                excl_rows(end+1,:) = {gName,rName,sub,taskTag,sprintf('BAD_NTRIALS_%d',height(Tuse))}; %#ok<SAGROW>
                continue;
            end

            [~,ord] = sort(Tuse.thisTrialN_raw,'ascend');
            Tuse = Tuse(ord,:);

            trialN = double(Tuse.thisTrialN_raw)';
            keepMask = (trialN >= CFG.DropTrials) & (trialN <= CFG.nTrials-1);

            keys_all = Tuse.(keyColName);
            keys_keep = keys_all(keepMask);

            if has_any_valid_key_response(keys_keep)
                keepSubs{end+1,1} = sub; %#ok<SAGROW>
            else
                excl_rows(end+1,:) = {gName,rName,sub,taskTag,'NO_KEYRESP_IN_USABLE_RUN'}; %#ok<SAGROW>
            end
        end

        FilteredMap.(gName).(rName) = keepSubs;
        fprintf('[FILTER] %s %s | original=%d | kept=%d | excluded=%d\n', ...
            gName, rName, numel(subjList), numel(keepSubs), numel(subjList)-numel(keepSubs));
    end
end

write_subjectlist_txt(FilteredMap, CFG.UpdatedSubjectListTxt);
fprintf('[WRITE] Updated subject list: %s\n', CFG.UpdatedSubjectListTxt);

if ~isempty(excl_rows)
    Texcl = cell2table(excl_rows, 'VariableNames', {'Group','Run','Subject','TaskTag','Reason'});
    writetable(Texcl, CFG.ExclusionLogCSV);
    fprintf('[WRITE] Exclusion log: %s\n', CFG.ExclusionLogCSV);
else
    fprintf('[INFO] No runs excluded.\n');
end

fprintf('\n====================================================\n');
fprintf('[STEP 2] Building RT-only conds_glm with TWO pmods using FILTERED subject list\n');
fprintf('====================================================\n');

qc_rows = {}; % Group Run Sub Model Nkeep N_nan N_abn RTcolFound N_ID_Missing N_Key_Missing N_Key_Valid

for gi = 1:size(PLAN,1)
    gName = PLAN{gi,1};
    runs  = PLAN{gi,2};

    for ri = 1:numel(runs)
        rName = runs{ri};
        taskTag = TaskTagMap.(gName).(rName);

        subjList = FilteredMap.(gName).(rName);
        behDir   = fullfile(CFG.BehaviorRoot, gName, rName);

        fprintf('\n----------------------------------------------------\n');
        fprintf('[CONDITIONS] %s %s | TaskTag=%s | nSub=%d\n', gName, rName, taskTag, numel(subjList));
        fprintf('----------------------------------------------------\n');

        outDir = fullfile(CFG.OutBase, sprintf('conds_glm_%s_%s_%s', gName, rName, CFG.ModelTag));
        if ~exist(outDir,'dir'); mkdir(outDir); end

        for si = 1:numel(subjList)
            sub = subjList{si};

            csv_file = fullfile(behDir, sprintf('%s_%s_%s_%s_trials.csv', ...
                gName, sub, lower(rName), taskTag));

            if ~exist(csv_file,'file')
                fprintf('[MISS] %s %s %s: %s\n', gName, rName, sub, csv_file);
                continue;
            end

            T = readtable(csv_file);

            if ~any(strcmp(T.Properties.VariableNames,'thisTrialN_raw'))
                error('Column thisTrialN_raw not found in %s', csv_file);
            end
            if ~any(strcmp(T.Properties.VariableNames,'Path_raw'))
                error('Column Path_raw not found in %s', csv_file);
            end

            [hasKeyCol, keyColName] = find_first_existing_var(T, CFG.KeyRespColCandidates);
            if ~hasKeyCol
                error('No key response column found in %s', csv_file);
            end

            use_mask = (T.thisTrialN_raw >= 0) & (T.thisTrialN_raw <= CFG.nTrials-1);
            Tuse = T(use_mask,:);
            if height(Tuse) ~= CFG.nTrials
                error('Expected %d trials (0..%d) but got %d in %s', ...
                    CFG.nTrials, CFG.nTrials-1, height(Tuse), csv_file);
            end

            [~,ord] = sort(Tuse.thisTrialN_raw,'ascend');
            Tuse = Tuse(ord,:);

            trialN = double(Tuse.thisTrialN_raw)';
            keepMask = (trialN >= CFG.DropTrials) & (trialN <= CFG.nTrials-1);
            trialN_keep = trialN(keepMask);
            nKeep = numel(trialN_keep);
            assert(nKeep==CFG.RetainedTrials,'Expected %d retained trials for %s %s %s.', ...
                CFG.RetainedTrials,gName,rName,sub);

            if nKeep <= 0
                error('After dropping first %d trial(s), no trials left for %s %s %s', ...
                    CFG.DropTrials, gName, rName, sub);
            end

            newIdx = trialN_keep - CFG.DropTrials;
            stim_onsets = newIdx * CFG.trialLen;
            rate_onsets = stim_onsets + CFG.stimDur;

            [rt_all, rtColFound] = get_numeric_column_or_nan(Tuse, CFG.RT_ColCandidates);
            rt_keep = rt_all(keepMask);

            valid = ~isnan(rt_keep) & (rt_keep > CFG.RT_Min) & (rt_keep < CFG.RT_Max);
            ratingDuration  = zeros(1,nKeep);
            ratingDuration(valid) = rt_keep(valid);
            fixationDuration = CFG.rateDur - ratingDuration;

            path_all = string(Tuse.Path_raw);
            path_keep = path_all(keepMask);

            img_id_keep = nan(1,nKeep);
            for kk = 1:nKeep
                img_id_keep(kk) = lookup_img_id(path_keep(kk), IMG);
            end
            n_id_missing = sum(isnan(img_id_keep));
            assert(n_id_missing==0,'Unmapped stimulus image in %s; inspect cfg.image_id_csv.',csv_file);

            order_keep = 1:nKeep;
            if CFG.OrderZScore
                order_param = zscore_safe(order_keep);
            else
                order_param = double(order_keep);
            end

            key_all = Tuse.(keyColName);
            key_keep = key_all(keepMask);
            [resp_bin, key_raw_str] = parse_key_response_bin(key_keep, CFG.RespMissingValue);
            n_key_missing = sum(resp_bin == CFG.RespMissingValue & ~strcmp(key_raw_str, '1') & ~strcmp(key_raw_str, '2'));
            n_key_valid   = sum(resp_bin == -1 | resp_bin == 1);

            names = {};
            onsets = {};
            durations = {};

            names{end+1}     = 'Stim';
            onsets{end+1}    = stim_onsets;
            durations{end+1} = repmat(CFG.stimDur, 1, nKeep);

            if CFG.IncludeRating
                names{end+1}     = 'Rating';
                onsets{end+1}    = rate_onsets;
                durations{end+1} = ratingDuration;
            end

            names{end+1}     = 'Fixation';
            onsets{end+1}    = rate_onsets + ratingDuration;
            durations{end+1} = fixationDuration;

            [names,onsets,durations] = sanitize_multicond(names,onsets,durations);

            if CFG.AddStimOrderPmod || CFG.AddStimRespPmod
                pmod = build_stim_two_pmods(names, CFG, order_param, resp_bin);
                orth = ones(1, numel(names));
                stimIdx = find(strcmp(names,'Stim'), 1);
                if ~isempty(stimIdx)
                    orth(stimIdx) = CFG.StimOrth;
                end
            end

            out_mat = fullfile(outDir, sprintf('Condition_GLM_%s_%s_%s_%s.mat', ...
                CFG.ModelTag, gName, sub, rName));

            if CFG.AddStimOrderPmod || CFG.AddStimRespPmod
                save(out_mat, 'names','onsets','durations','pmod','orth');
            else
                save(out_mat, 'names','onsets','durations');
            end

            if CFG.WriteStimMetaCSV
                Tmeta = table();
                Tmeta.BetaIndex      = (1:nKeep)';
                Tmeta.OriginalTrial  = trialN_keep(:)+1; % Human-readable 1-based trial.
                Tmeta.ImageName      = string(arrayfun(@(x)char(basename_key(x)),path_keep(:),'UniformOutput',false));
                Tmeta.thisTrialN_raw = trialN_keep(:);
                Tmeta.NewIdx0        = newIdx(:);
                Tmeta.StimOnset      = stim_onsets(:);
                Tmeta.Order          = order_keep(:);
                Tmeta.OrderParam     = order_param(:);
                Tmeta.KeyRaw         = string(key_raw_str(:));
                Tmeta.RespBin        = resp_bin(:);
                Tmeta.Path_raw       = path_keep(:);
                Tmeta.Img_ID         = img_id_keep(:);

                out_meta = fullfile(outDir, sprintf('StimMeta_%s_%s_%s_%s.csv', ...
                    CFG.ModelTag, gName, sub, rName));
                writetable(Tmeta, out_meta);
            end

            n_nan = sum(isnan(rt_keep));
            n_abn = sum(~isnan(rt_keep) & (rt_keep<=CFG.RT_Min | rt_keep>=CFG.RT_Max));

            qc_rows(end+1,:) = {gName,rName,sub,CFG.ModelTag,nKeep,n_nan,n_abn,rtColFound,n_id_missing,n_key_missing,n_key_valid}; %#ok<SAGROW>

            condEnd = max([ ...
                max(stim_onsets + CFG.stimDur), ...
                max((rate_onsets + ratingDuration) + fixationDuration) ...
            ]);
            scanLen = nKeep * CFG.trialLen;
            fprintf('[OK] %s %s %s | nKeep=%d | validKey=%d | missKey=%d | missID=%d | scan=%.1fs | condEnd=%.3fs | RTcol=%s\n', ...
                gName,rName,sub,nKeep,n_key_valid,n_key_missing,n_id_missing,scanLen,condEnd,rtColFound);
        end
    end
end

if ~isempty(qc_rows)
    Tqc = cell2table(qc_rows, 'VariableNames', ...
        {'Group','Run','Subject','Model','Nkeep','N_RT_NaN','N_RT_Abnormal','RT_Column','N_ImgID_Missing','N_Key_Missing','N_Key_Valid'});
    writetable(Tqc, CFG.QCSummaryCSV);
    fprintf('\n[WRITE] %s\n', CFG.QCSummaryCSV);
end

fprintf('\n[ALL DONE] RT-only conds_glm build with TWO pmods finished.\n');

end

function M = tasktag_map()
    M = struct();
    M.Group1 = struct('Run1','sa','Run2','li','Run3','bea','Run4','img');
    M.Group2 = struct('Run1','li','Run2','sa','Run3','img','Run4','bea');
    M.Group3 = struct('Run1','bea','Run2','img','Run3','sa','Run4','li');
    M.Group4 = struct('Run1','img','Run2','bea','Run3','li','Run4','sa');
end

function [tf, varname] = find_first_existing_var(T, candidates)
    tf = false;
    varname = '';
    vn = T.Properties.VariableNames;
    for i = 1:numel(candidates)
        c = candidates{i};
        c2 = strrep(c,'.','_');
        tryList = {c, c2};
        for j = 1:numel(tryList)
            cc = tryList{j};
            if any(strcmp(vn, cc))
                tf = true;
                varname = cc;
                return;
            end
        end
    end
end

function tf = has_any_valid_key_response(v)
    tf = false;
    if isempty(v), return; end

    if isnumeric(v)
        tf = any(~isnan(v(:)));
        return;
    end

    if iscell(v)
        for i = 1:numel(v)
            x = v{i};
            if isnumeric(x)
                if ~isempty(x) && any(~isnan(x))
                    tf = true; return;
                end
            else
                s = strtrim(char(string(x)));
                if ~isempty(s) && ~strcmpi(s,'nan') && ~strcmpi(s,'none') && ~strcmpi(s,'null') ...
                        && ~strcmpi(s,'[]')
                    tf = true; return;
                end
            end
        end
        return;
    end

    if isstring(v) || ischar(v)
        vs = string(v);
        vs = strtrim(vs);
        tf = any(vs ~= "" & lower(vs) ~= "nan" & lower(vs) ~= "none" & lower(vs) ~= "null" & vs ~= "[]");
        return;
    end
end

function [resp_bin, raw_str] = parse_key_response_bin(v, missingValue)

    n = numel(v);
    resp_bin = repmat(missingValue, 1, n);
    raw_str  = repmat({''}, 1, n);

    for i = 1:n
        x = [];
        s = '';

        if iscell(v)
            xi = v{i};
        else
            xi = v(i);
        end

        if isnumeric(xi)
            if isempty(xi) || any(isnan(xi))
                x = NaN;
                s = '';
            else
                x = double(xi(1));
                s = strtrim(num2str(x));
            end
        else
            s = strtrim(char(string(xi)));
            if isempty(s) || strcmpi(s,'nan') || strcmpi(s,'none') || strcmpi(s,'null') || strcmpi(s,'[]')
                x = NaN;
                s = '';
            else
                tok = regexp(s, '(-?\d+)', 'tokens', 'once');
                if isempty(tok)
                    x = NaN;
                else
                    x = str2double(tok{1});
                end
            end
        end

        raw_str{i} = s;
        if ~isnan(x)
            if x == 1
                resp_bin(i) = -1;
            elseif x == 2
                resp_bin(i) = 1;
            else
                resp_bin(i) = missingValue;
            end
        end
    end
end

function [x, usedCol] = get_numeric_column_or_nan(T, candidates)
    usedCol = '';
    x = nan(height(T),1);
    vars = T.Properties.VariableNames;

    for k = 1:numel(candidates)
        c = candidates{k};
        c2 = strrep(c,'.','_');
        tryList = {c, c2};

        for j = 1:numel(tryList)
            cc = tryList{j};
            if any(strcmp(vars, cc))
                v = T.(cc);
                if iscell(v)
                    v = cellfun(@str2double, v);
                end
                x = double(v);
                usedCol = cc;
                return;
            end
        end
    end
    usedCol = '(RT col missing -> all NaN)';
end

function [names,onsets,durations] = sanitize_multicond(names,onsets,durations)
% Remove nonfinite event entries while keeping onsets and durations paired.
    n = numel(names);
    keep = true(1,n);

    for i = 1:n
        if isempty(onsets{i})
            keep(i) = false;
            continue;
        end
        onsets{i} = double(onsets{i}(:))';
        durations{i} = double(durations{i}(:))';
        if numel(durations{i}) == 1 && numel(onsets{i}) > 1
            durations{i} = repmat(durations{i}, 1, numel(onsets{i}));
        end
        if numel(durations{i}) ~= numel(onsets{i})
            error('Condition "%s": duration length (%d) != onset length (%d)', ...
                names{i}, numel(durations{i}), numel(onsets{i}));
        end
    end

    if any(~keep)
        fprintf('[WARN] Dropping empty-onset conditions: %s\n', strjoin(names(~keep),', '));
        names = names(keep);
        onsets = onsets(keep);
        durations = durations(keep);
    end
end

function IMG = load_imageid_maps(csvPath, fileColCands, idColCands)
    if ~exist(csvPath,'file')
        error('imageID.csv not found: %s', csvPath);
    end
    T = readtable(csvPath, 'TextType','string');

    fileCol = pick_col(T, fileColCands);
    idCol   = pick_col(T, idColCands);

    file_raw = string(T.(fileCol));
    id_raw   = T.(idCol);

    id_num = double(id_raw);
    assert(numel(unique(id_num))==numel(id_num) && all(isfinite(id_num)), ...
        'Img_ID values must be unique finite stimulus indices.');
    ok = ~ismissing(file_raw) & ~isnan(id_num);

    file_raw = file_raw(ok);
    id_num   = id_num(ok);

    assert(numel(unique(lower(file_raw)))==numel(file_raw),'Duplicate File_Name in image ID table.');
    key_full  = normalize_key(file_raw);
    key_base  = basename_key(file_raw);
    key_noext = basename_noext_key(file_raw);

    IMG = struct();
    IMG.fileCol = fileCol;
    IMG.idCol   = idCol;
    IMG.N       = numel(id_num);

    IMG.map_full  = make_map_unique(key_full,  id_num, 'FULL');
    IMG.map_base  = make_map_unique(key_base,  id_num, 'BASE');
    IMG.map_noext = make_map_unique(key_noext, id_num, 'NOEXT');
end

function col = pick_col(T, candidates)
    vn = string(T.Properties.VariableNames);
    col = '';
    for i = 1:numel(candidates)
        c = string(candidates{i});
        idx = find(strcmpi(vn, c), 1);
        if ~isempty(idx)
            col = char(vn(idx));
            return;
        end
    end
    error('Cannot find required column. Candidates: %s', strjoin(string(candidates), ', '));
end

function mp = make_map_unique(keys, vals, tag)
    keys = string(keys);
    vals = double(vals);

    [uk, ia] = unique(keys, 'stable');
    if numel(uk) < numel(keys)
        fprintf('[IMG][WARN] %s keys duplicated: %d -> %d (keep first)\n', tag, numel(keys), numel(uk));
    end
    uk = uk(:);
    uv = vals(ia);

    mp = containers.Map('KeyType','char','ValueType','double');
    for i = 1:numel(uk)
        mp(char(uk(i))) = uv(i);
    end
end

function id = lookup_img_id(pathStr, IMG)
    id = nan;
    if ismissing(pathStr) || strlength(pathStr)==0
        return;
    end

    k_full  = char(normalize_key(pathStr));
    k_base  = char(basename_key(pathStr));
    k_noext = char(basename_noext_key(pathStr));

    if isKey(IMG.map_full, k_full)
        id = IMG.map_full(k_full); return;
    end
    if isKey(IMG.map_base, k_base)
        id = IMG.map_base(k_base); return;
    end
    if isKey(IMG.map_noext, k_noext)
        id = IMG.map_noext(k_noext); return;
    end
end

function k = normalize_key(p)
    p = string(p);
    p = strtrim(p);
    p = replace(p, "\\", "/");
    p = lower(p);
    k = p;
end

function k = basename_key(p)
    p = normalize_key(p);
    parts = split(p, "/");
    b = parts(end);
    k = b;
end

function k = basename_noext_key(p)
    b = basename_key(p);
    [~, name, ~] = fileparts(char(b));
    k = lower(string(name));
end

function z = zscore_safe(x)
    x = double(x(:));
    mu = mean(x);
    sd = std(x);
    if sd < eps
        z = zeros(size(x));
    else
        z = (x - mu) ./ sd;
    end
    z = z(:)';
end

function pmod = build_stim_two_pmods(names, CFG, order_param, resp_bin)
% Trial order and signed response modulate Stim independently (orth=0).
    ncond = numel(names);
    pmod = repmat(struct('name',{{}}, 'param',{{}}, 'poly',{{}}), 1, ncond);

    stimIdx = find(strcmp(names,'Stim'), 1);
    if isempty(stimIdx)
        error('Cannot find condition "Stim" in names.');
    end

    nm = {};
    pm = {};
    pl = {};

    if CFG.AddStimOrderPmod
        nm{end+1} = 'Order'; %#ok<AGROW>
        pm{end+1} = double(order_param(:)'); %#ok<AGROW>
        pl{end+1} = double(CFG.OrderPoly); %#ok<AGROW>
    end

    if CFG.AddStimRespPmod
        nm{end+1} = 'RespBin'; %#ok<AGROW>
        pm{end+1} = double(resp_bin(:)'); %#ok<AGROW>
        pl{end+1} = double(CFG.RespPoly); %#ok<AGROW>
    end

    pmod(stimIdx).name  = nm;
    pmod(stimIdx).param = pm;
    pmod(stimIdx).poly  = pl;
end

function SubjectMap = parse_subjectlist_txt(txtPath)
    if ~exist(txtPath,'file')
        error('SubjectList.txt not found: %s', txtPath);
    end

    raw = fileread(txtPath);
    raw = strrep(raw, sprintf('\r\n'), sprintf('\n'));
    raw = strrep(raw, sprintf('\r'),   sprintf('\n'));
    lines = regexp(raw, '\n', 'split');
    lines = lines(:);

    clean = cell(size(lines));
    for i = 1:numel(lines)
        s = char(string(lines{i}));
        s = strrep(s, char(65279), '');
        s = strrep(s, char(160), ' ');
        s = strrep(s, sprintf('\t'), ' ');
        s = strtrim(s);
        clean{i} = s;
    end

    SubjectMap = struct();
    groups = {'Group1','Group2','Group3','Group4'};
    runs   = {'Run1','Run2','Run3','Run4'};

    for g = 1:numel(groups)
        gName = groups{g};
        gIdx = find(strcmpi(clean, gName), 1, 'first');
        if isempty(gIdx)
            SubjectMap.(gName) = struct('Run1',{{}},'Run2',{{}},'Run3',{{}},'Run4',{{}});
            continue;
        end

        endIdx = numel(clean)+1;
        for i = gIdx+1:numel(clean)
            if ~isempty(regexp(clean{i}, '^Group\d+$', 'once'))
                endIdx = i; break;
            end
        end
        block = clean(gIdx+1:endIdx-1);

        SubjectMap.(gName) = struct();
        for r = 1:numel(runs)
            rName = runs{r};
            runStart = [];
            for i = 1:numel(block)
                if ~isempty(regexp(block{i}, ['^' regexptranslate('escape', rName) '\s*='], 'once'))
                    runStart = i; break;
                end
            end
            if isempty(runStart)
                SubjectMap.(gName).(rName) = {};
                continue;
            end

            buf = block{runStart};
            j = runStart;
            while isempty(strfind(buf, '};')) && j < numel(block)
                j = j + 1;
                buf = [buf ' ' block{j}]; %#ok<AGROW>
            end

            tok = regexp(buf, ['^' regexptranslate('escape', rName) '\s*=\s*\{([^}]*)\}\s*;'], 'tokens', 'once');
            if isempty(tok) || isempty(tok{1})
                SubjectMap.(gName).(rName) = {};
                continue;
            end
            inside = tok{1};

            idTok = regexp(inside, '''(\d+)''', 'tokens');
            subjList = cell(numel(idTok), 1);
            for k = 1:numel(idTok)
                sid = idTok{k}{1};
                subjList{k} = sprintf('%03d', str2double(sid));
            end
            SubjectMap.(gName).(rName) = subjList;
        end
    end
end

function write_subjectlist_txt(SubjectMap, outPath)
    fid = fopen(outPath, 'w');
    if fid == -1
        error('Cannot open output subject list: %s', outPath);
    end

    groups = {'Group1','Group2','Group3','Group4'};
    runs   = {'Run1','Run2','Run3','Run4'};

    for g = 1:numel(groups)
        gName = groups{g};
        fprintf(fid, '%s\n', gName);

        for r = 1:numel(runs)
            rName = runs{r};
            subjList = SubjectMap.(gName).(rName);

            fprintf(fid, '%s = {', rName);
            for i = 1:numel(subjList)
                fprintf(fid, '''%s''', subjList{i});
                if i < numel(subjList)
                    fprintf(fid, ',');
                end
            end
            fprintf(fid, '};\n');
        end
        fprintf(fid, '\n');
    end

    fclose(fid);
end
