function run_beta_series_connectivity(cfg)
% RUN_BETA_SERIES_CONNECTIVITY Compute 15-ROI beta-series connectivity.
% Atlas labels are sampled in the beta-image grid by nearest neighbour.
% Trial-wise voxel means form a 99-by-15 beta series. Pearson correlations
% and Fisher-z transforms provide 105 unique connections per subject/task.
% Outputs include individual matrices, a complete-case edge table and Zall
% indexed as ROI-by-ROI-by-subject-by-task.

if nargin < 1, cfg = study_config(); end
fields = {'spm_dir','lss_root','subject_list','merged_atlas_dir','hcpex_atlas','hcpex_labels','connectivity_root'};
for k = 1:numel(fields)
    assert(isfield(cfg,fields{k}) && ~isempty(cfg.(fields{k})), 'Missing cfg.%s',fields{k});
end
assert(isfolder(cfg.lss_root),'Missing LSS root: %s',cfg.lss_root);
assert(isfile(cfg.subject_list),'Missing subject list: %s',cfg.subject_list);
addpath(cfg.spm_dir);
spm('Defaults','fMRI');

opt = struct('ExpectedTrialsPerRun',99,'MinimumTrials',3,'Overwrite',false);
if isfield(cfg,'n_trials'), opt.ExpectedTrialsPerRun = cfg.n_trials; end
if isfield(cfg,'connectivity_options')
    keys = fieldnames(cfg.connectivity_options);
    for k = 1:numel(keys)
        assert(isfield(opt,keys{k}),'Unknown connectivity_options key: %s',keys{k});
        opt.(keys{k}) = cfg.connectivity_options.(keys{k});
    end
end
tasks = {'Safety','Vitality','Beauty','Imageability'};
PLAN = cell(4,2);
for g = 1:4, PLAN(g,:) = {sprintf('Group%d',g),{'Run1','Run2','Run3','Run4'}}; end
if isfield(cfg,'plan'), PLAN = cfg.plan; end
root = cfg.connectivity_root;
if ~isfolder(root), mkdir(root); end

records = struct('Group',{},'Run',{},'Subject',{},'Task',{},'BetaDir',{});
for g = 1:size(PLAN,1)
    group = PLAN{g,1}; runs = PLAN{g,2};
    for r = 1:numel(runs)
        run = runs{r}; task = task_from_group_run(group,run);
        subjects = load_subjectlist_from_txt_exact(group,run,cfg.subject_list);
        for s = 1:numel(subjects)
            subject = subjects{s};
            records(end+1) = struct('Group',group,'Run',run,'Subject',subject, ...
                'Task',task,'BetaDir',fullfile(cfg.lss_root,group,run, ...
                ['UP_' subject],['LSS_StimBeta_' run])); %#ok<AGROW>
        end
    end
end
assert(~isempty(records),'No records found in subject list.');
assert(numel(unique(string({records.Group})+"|"+string({records.Subject})+"|"+ ...
    string({records.Task})))==numel(records),'Duplicate subject x task records.');

firstBeta = '';
for i = 1:numel(records)
    dd = dir(fullfile(records(i).BetaDir,'beta_stim_trial_*.nii'));
    if ~isempty(dd), firstBeta = fullfile(dd(1).folder,dd(1).name); break; end
end
assert(~isempty(firstBeta),'No beta_stim_trial_*.nii files found in valid runs.');
Atlas = prepare_atlas(cfg,spm_vol(firstBeta));
save(fullfile(root,'PreparedAtlas_HCPex15.mat'),'Atlas','-v7.3');
writetable(table(Atlas.roi_ids,string(Atlas.roi_labels),cellfun(@numel,Atlas.atlas_idx), ...
    'VariableNames',{'ROI_ID','ROI','NVoxels'}),fullfile(root,'PreparedAtlas_HCPex15_Lookup.csv'));

edgeRows = {}; qcRows = {}; results = cell(numel(records),1);
[edgeI,edgeJ] = find(triu(true(15),1));
for i = 1:numel(records)
    rec = records(i);
    dd = dir(fullfile(rec.BetaDir,'beta_stim_trial_*.nii'));
    if isempty(dd)
        qcRows(end+1,:) = {rec.Group,rec.Subject,rec.Run,rec.Task,0, ...
            opt.ExpectedTrialsPerRun,false,'Missing LSS beta files'}; %#ok<AGROW>
        continue;
    end
    trial_numbers = zeros(numel(dd),1);
    for t = 1:numel(dd)
        tok = regexp(dd(t).name,'beta_stim_trial_(\d+)\.nii$','tokens','once');
        assert(~isempty(tok),'Unrecognized beta file: %s',dd(t).name);
        trial_numbers(t) = str2double(tok{1});
    end
    [trial_numbers,order] = sort(trial_numbers); dd = dd(order);
    assert(numel(unique(trial_numbers))==numel(dd),'Duplicate beta indices: %s',rec.BetaDir);
    nTrial = numel(dd);
    assert(nTrial==opt.ExpectedTrialsPerRun,'Expected %d beta images in %s.', ...
        opt.ExpectedTrialsPerRun,rec.BetaDir);
    assert(nTrial>=opt.MinimumTrials,'Too few LSS trials in %s',rec.BetaDir);
    seriesDir = fullfile(root,'Step1_ROIBetaSeries',rec.Task);
    bscDir = fullfile(root,'Step2_BSCMatrices',rec.Task);
    if ~isfolder(seriesDir), mkdir(seriesDir); end
    if ~isfolder(bscDir), mkdir(bscDir); end
    basename = sprintf('%s_UP_%s_%s_HCPex15',rec.Group,rec.Subject,rec.Run);
    seriesFile = fullfile(seriesDir,[basename '_betaSeries.mat']);
    beta_files = string(fullfile({dd.folder},{dd.name}))';
    if isfile(seriesFile) && ~opt.Overwrite
        old = load(seriesFile,'X','trial_numbers','roi_ids','beta_files', ...
            'source_sizes','source_dates','atlas_signature');
        assert(all(isfield(old,{'X','trial_numbers','roi_ids','beta_files', ...
            'source_sizes','source_dates','atlas_signature'})) && ...
            isequal(size(old.X),[nTrial 15]) && isequal(old.beta_files,beta_files) && ...
            isequal(old.trial_numbers,trial_numbers) && isequal(old.roi_ids,Atlas.roi_ids) && ...
            isequal(old.source_sizes,[dd.bytes]') && isequal(old.source_dates,[dd.datenum]') && ...
            isfield(old,'atlas_signature') && isequal(old.atlas_signature,Atlas.signature), ...
            'Cached inputs changed. Set connectivity_options.Overwrite=true: %s',seriesFile);
        X = old.X;
    else
        X = nan(nTrial,15);
        for t = 1:nTrial
            V = spm_vol(char(beta_files(t)));
            assert(isequal(V.dim,Atlas.ref_dim) && max(abs(V.mat(:)-Atlas.ref_mat(:)))<1e-6, ...
                'Beta and atlas grid mismatch: %s',beta_files(t));
            Y = spm_read_vols(V);
            for r = 1:15
                values = Y(Atlas.atlas_idx{r}); values = values(isfinite(values));
                if ~isempty(values), X(t,r) = mean(values); end
            end
        end
        roi_ids = Atlas.roi_ids; roi_labels = Atlas.roi_labels; %#ok<NASGU>
        source_sizes = [dd.bytes]'; source_dates = [dd.datenum]'; %#ok<NASGU>
        atlas_signature = Atlas.signature; %#ok<NASGU>
        save(seriesFile,'X','trial_numbers','roi_ids','roi_labels','rec','beta_files', ...
            'source_sizes','source_dates','atlas_signature','-v7.3');
    end
    [R,Z,bad_roi] = beta_series_matrix(X);
    roi_ids = Atlas.roi_ids; roi_labels = Atlas.roi_labels; %#ok<NASGU>
    save(fullfile(bscDir,[basename '_BSC.mat']),'R','Z','bad_roi', ...
        'rec','trial_numbers','roi_ids','roi_labels','-v7.3');
    results{i} = Z;
    subjectKey = string(rec.Group)+"_"+string(rec.Subject);
    for e = 1:numel(edgeI)
        a = edgeI(e); b = edgeJ(e);
        edgeRows(end+1,:) = {char(subjectKey),rec.Subject,rec.Group,rec.Run,rec.Task,e, ...
            Atlas.roi_ids(a),Atlas.roi_labels{a},Atlas.roi_ids(b),Atlas.roi_labels{b}, ...
            R(a,b),Z(a,b),nTrial}; %#ok<AGROW>
    end
    qcRows(end+1,:) = {rec.Group,rec.Subject,rec.Run,rec.Task,nTrial, ...
        opt.ExpectedTrialsPerRun,nTrial==opt.ExpectedTrialsPerRun, ...
        sprintf('%d invalid ROIs',sum(bad_roi))}; %#ok<AGROW>
end
assert(~isempty(edgeRows),'No usable beta-series runs found.');
T = cell2table(edgeRows,'VariableNames',{'Subject','OriginalSubject','Group','Run','Task', ...
    'Edge','ROI1_ID','ROI1','ROI2_ID','ROI2','R','Z','NTrials'});
QC = cell2table(qcRows,'VariableNames',{'Group','Subject','Run','Task','NTrials', ...
    'ExpectedTrials','MatchesExpected','Status'});
writetable(QC,fullfile(root,'BSC_Run_QC.csv'));

subject_keys = unique(string(T.Subject),'stable'); keep = false(numel(subject_keys),1);
for s = 1:numel(subject_keys)
    keep(s) = numel(unique(string(T.Task(string(T.Subject)==subject_keys(s)))))==4;
end
subjects = subject_keys(keep);
assert(~isempty(subjects),'No subjects have all four tasks. Inspect BSC_Run_QC.csv.');
included = ismember(string(T.Subject),subjects);
writetable(T,fullfile(root,'BSC_EdgeLong_AllAvailable.csv'));
T = T(included,:);
writetable(T,fullfile(root,'BSC_EdgeLong.csv'));
Zall = nan(15,15,numel(subjects),4);
for s = 1:numel(subjects)
    for t = 1:4
        idx = find((string({records.Group})+"_"+string({records.Subject})==subjects(s)) & ...
            string({records.Task})==tasks{t});
        assert(isscalar(idx) && ~isempty(results{idx}),'Ambiguous complete-case run.');
        Zall(:,:,s,t) = results{idx};
    end
end
roi_ids = Atlas.roi_ids; roi_labels = Atlas.roi_labels; %#ok<NASGU>
save(fullfile(root,'BSC_all_tasks_stack_HCPex15.mat'),'Zall','subjects','tasks', ...
    'roi_ids','roi_labels','-v7.3');
fprintf('[BSC] %d complete-case subjects; %d connections; QC: %s\n', ...
    numel(subjects),numel(edgeI),fullfile(root,'BSC_Run_QC.csv'));
end

function Atlas = prepare_atlas(cfg,Vref)
atlasFile = cfg.hcpex_atlas;
lookupFile = cfg.hcpex_labels;
assert(isfile(atlasFile) && isfile(lookupFile),'Build the merged atlas first.');
LUT = readtable(lookupFile,'Encoding','UTF-8','VariableNamingRule','preserve');
[ids,first] = unique(double(LUT.Merged_ROI_ID),'stable');
labels = cellstr(string(LUT.Merged_ROI_Abbrev(first)));
normalized = lower(regexprep(string(labels),'[^a-zA-Z0-9]',''));
keep = ~(contains(normalized,'amyg') | contains(normalized,'entorhinal') | normalized=="ec");
ids = ids(keep); labels = labels(keep);
assert(numel(ids)==15,'Expected 15 retained ROIs after excluding amygdala and entorhinal.');
Vsrc = spm_vol(atlasFile); A = zeros(Vref.dim);
[x,y] = ndgrid(1:Vref.dim(1),1:Vref.dim(2));
for z = 1:Vref.dim(3)
    reference = [x(:)';y(:)';repmat(z,1,numel(x));ones(1,numel(x))];
    source = Vsrc.mat \ (Vref.mat*reference);
    values = spm_sample_vol(Vsrc,source(1,:),source(2,:),source(3,:),0);
    A(:,:,z) = reshape(round(values),size(x));
end
indices = cell(15,1);
for r = 1:15
    indices{r} = find(A==ids(r));
    assert(~isempty(indices{r}),'ROI %s is empty in beta space.',labels{r});
end
Atlas = struct('roi_ids',ids(:),'roi_labels',{labels(:)},'atlas_idx',{indices}, ...
    'ref_dim',Vref.dim,'ref_mat',Vref.mat,'source_atlas',atlasFile);
atlasInfo = dir(atlasFile); lookupInfo = dir(lookupFile);
Atlas.signature = struct('atlas_file',atlasFile,'atlas_size',atlasInfo.bytes, ...
    'atlas_modified',atlasInfo.datenum,'lookup_file',lookupFile, ...
    'lookup_size',lookupInfo.bytes,'lookup_modified',lookupInfo.datenum, ...
    'reference_dim',Vref.dim,'reference_affine',Vref.mat);
end

function task = task_from_group_run(group,run)
map = {'Safety','Vitality','Beauty','Imageability'; ...
    'Vitality','Safety','Imageability','Beauty'; ...
    'Beauty','Imageability','Safety','Vitality'; ...
    'Imageability','Beauty','Vitality','Safety'};
g = sscanf(group,'Group%d'); r = sscanf(run,'Run%d');
assert(isscalar(g) && isscalar(r) && ismember(g,1:4) && ismember(r,1:4), ...
    'Unrecognized group/run: %s %s',group,run);
task = map{g,r};
end

function subjects = load_subjectlist_from_txt_exact(group,run,path)
raw = regexprep(fileread(path),'\r\n?','\n');
raw = strrep(raw,char(65279),'');
lines = regexp(raw,'\n','split'); lines = strtrim(lines);
start = find(strcmp(lines,group),1);
assert(~isempty(start),'Missing %s in subject list.',group);
stop = numel(lines)+1;
for i = start+1:numel(lines)
    if ~isempty(regexp(lines{i},'^Group\d+$','once')), stop = i; break; end
end
block = strjoin(lines(start+1:stop-1),newline);
token = regexp(block,['(?m)^' regexptranslate('escape',run) '\s*=\s*\{([^}]*)\}\s*;'], ...
    'tokens','once');
assert(~isempty(token),'Missing %s %s subject list.',group,run);
ids = regexp(token{1},'''(\d+)''','tokens');
subjects = cellfun(@(x)sprintf('%03d',str2double(x{1})),ids,'UniformOutput',false);
assert(numel(unique(subjects))==numel(subjects),'Duplicate subjects in %s %s.',group,run);
end
