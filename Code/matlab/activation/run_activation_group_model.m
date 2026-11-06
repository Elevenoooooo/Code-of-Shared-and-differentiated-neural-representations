function run_activation_group_model(cfg)
% RUN_ACTIVATION_GROUP_MODEL Estimate task activation and task differences.
% Positive activation: a separate pooled one-sample T model for each task.
% Task differences: repeated-measures Subject + GroupTaskCell parameterization.
% Equal-group marginal contrasts compare Safety, Vitality, Beauty, Imageability.
% One signed A-B T map is retained per pair, with two-sided voxelwise FDR.
% ROI labels are sampled in the first-level functional grid by nearest neighbour.
% Outputs include unthresholded statistics, p/q maps, masks and ROI summaries.
if nargin < 1, cfg = study_config(); end
required = {'spm_dir','first_level_root','activation_root','merged_atlas_dir','hcpex_atlas','hcpex_labels'};
for k = 1:numel(required)
    assert(isfield(cfg,required{k}) && ~isempty(cfg.(required{k})), 'Missing cfg.%s', required{k});
end
addpath(cfg.spm_dir);

spm('Defaults','fMRI');
spm_jobman('initcfg');

CFG.FirstLevelRoot = cfg.first_level_root;

CFG.OutRoot = cfg.activation_root;

CFG.AtlasDir = cfg.merged_atlas_dir;
CFG.AtlasNii = cfg.hcpex_atlas;
CFG.LookupCsv = cfg.hcpex_labels;
CFG.RoiCountCsv = fullfile(CFG.AtlasDir, 'HCPex_MergedROI_UrbanPerception_Refined17_roi_voxelcount.csv');

CFG.UnionMaskNii = fullfile(CFG.AtlasDir, 'HCPex_MergedROI_UrbanPerception_Refined17_unionmask_2mm.nii');

CFG.Groups = {'Group1','Group2','Group3','Group4'};
CFG.Runs   = {'Run1','Run2','Run3','Run4'};
CFG.Tasks  = {'Safety','Vitality','Beauty','Imageability'};

CFG.ContrastSets = { ...
    struct('Name','Stim_gt_Fixation', 'FirstLevelConName','Stim>Fixation'), ...
    };

CFG.PairList = { ...
    'Safety',   'Vitality'; ...
    'Safety',   'Beauty'; ...
    'Safety',   'Imageability'; ...
    'Vitality', 'Beauty'; ...
    'Vitality', 'Imageability'; ...
    'Beauty',   'Imageability' ...
    };

CFG.EqualGroupWeight = true;
CFG.DebugContrastColumns = false;

CFG.AddInteractionF = false;

CFG.WriteBidirectionalPairwise = false;
CFG.PairwiseT_TwoSidedFDR = true;

CFG.AddPairwiseEqualityF = false;
CFG.PairwiseEqualityF_Suffix = 'neq_F';


CFG.RequireCompleteSubject = true;

CFG.CleanAndResliceConImages = true;
CFG.ResliceInterp = 1;

CFG.RecreateContrasts = true;

CFG.FDR_q = 0.05;
CFG.WriteVoxelCSV = true;
CFG.WriteROISummaryCSV = true;
CFG.WritePmapQmapSigmask = true;
CFG.WriteStatMap = true;
CFG.WriteThresholdedStatMap = true;

CFG.FDR_IncludeTaskOneSample = true;   % include separate task one-sample maps in STEP4 FDR outputs
CFG.FDR_IncludeInteractionF  = false;

CFG.FDR_CorrectionScope      = 'ROI';
CFG.FDR_OneSample_CorrectionScope = 'Union';
CFG.WriteFDRByROIFolders     = true;   % only TaskMainEffect / Pairwise / PairwiseEqualityF are split into ByROI folders

CFG.RunMode = 'full';

CFG.OverwriteModels = false;

CFG.ContrastSetsToRun = {'Stim_gt_Fixation'};

CFG.SkipExistingAtlasConImages = true;

CFG.AutoRebuildManifestIfInvalid = true;

CFG.UseExplicitROIMaskInSPMModel = false;

CFG.RunPreflightVoxelQC = true;
CFG.MinNonzeroVoxInROIMask = 10;

CFG.MaxVoxelRowsPerContrast = Inf;

if isfield(cfg,'activation_group_options')
    user = cfg.activation_group_options; names = fieldnames(user);
    for k = 1:numel(names)
        assert(isfield(CFG,names{k}), 'Unknown activation_group_options option: %s',names{k});
        CFG.(names{k}) = user.(names{k});
    end
end
switch lower(CFG.RunMode)
    case 'full'
        CFG.RunStep1_BuildManifest       = true;
        CFG.RunStep2_CleanReslice        = true;
        CFG.RunStep3A_TaskOneSample      = true;
        CFG.RunStep3_EstimateSecondLevel = true;
        CFG.RunStep3_AddContrasts        = true;
        CFG.RunStep4_FDROutputs          = true;

    case 'models_only'
        CFG.RunStep1_BuildManifest       = false;
        CFG.RunStep2_CleanReslice        = false;
        CFG.RunStep3A_TaskOneSample      = true;
        CFG.RunStep3_EstimateSecondLevel = true;
        CFG.RunStep3_AddContrasts        = true;
        CFG.RunStep4_FDROutputs          = true;

    case 'outputs_only'
        CFG.RunStep1_BuildManifest       = false;
        CFG.RunStep2_CleanReslice        = false;
        CFG.RunStep3A_TaskOneSample      = false;
        CFG.RunStep3_EstimateSecondLevel = false;
        CFG.RunStep3_AddContrasts        = false;
        CFG.RunStep4_FDROutputs          = true;

    otherwise
        error('Unknown CFG.RunMode: %s. Use full, models_only, or outputs_only.', CFG.RunMode);
end

CFG.OverwriteTaskOneSample       = CFG.OverwriteModels;
CFG.OverwriteSecondLevel         = CFG.OverwriteModels;

assert(exist(CFG.FirstLevelRoot,'dir') == 7, 'Missing first-level root: %s', CFG.FirstLevelRoot);
assert(exist(CFG.AtlasNii,'file') == 2, 'Missing atlas NIfTI: %s', CFG.AtlasNii);
assert(exist(CFG.LookupCsv,'file') == 2, 'Missing lookup CSV: %s', CFG.LookupCsv);

CFG = prepare_functional_roi_atlas(CFG);

mkdir_if_needed(CFG.OutRoot);
mkdir_if_needed(fullfile(CFG.OutRoot, 'Manifest'));
mkdir_if_needed(fullfile(CFG.OutRoot, 'ConAtlasSpace'));
mkdir_if_needed(fullfile(CFG.OutRoot, 'SecondLevel_TaskOneSample'));
mkdir_if_needed(fullfile(CFG.OutRoot, 'SecondLevel_FullModel'));
mkdir_if_needed(fullfile(CFG.OutRoot, 'FDR_VoxelTables'));
mkdir_if_needed(fullfile(CFG.OutRoot, 'FDR_ROISummary'));
mkdir_if_needed(fullfile(CFG.OutRoot, 'FDR_Maps'));

CFG.ContrastSets = filter_contrast_sets(CFG.ContrastSets, CFG.ContrastSetsToRun);

ensure_union_mask(CFG);
LUT = load_roi_lookup(CFG);

manifestPath = fullfile(CFG.OutRoot, 'Manifest', 'FirstLevel_ContrastManifest.csv');
manifestAtlasPath = fullfile(CFG.OutRoot, 'Manifest', 'FirstLevel_ContrastManifest_atlas.csv');

if CFG.RunStep1_BuildManifest
    fprintf('\n[STEP1] Build manifest from first-level results...\n');
    Tall = build_manifest(CFG);
    Tall = standardize_manifest_vars(Tall, 'STEP1 build_manifest');
    write_table_utf8bom(Tall, manifestPath);
    fprintf('[WRITE] %s\n', manifestPath);
else
    fprintf('\n[STEP1-SKIP] Read existing manifest...\n');
    try
        assert(exist(manifestPath,'file') == 2, 'Missing manifest CSV: %s', manifestPath);
        Tall = safe_readtable(manifestPath);
        Tall = standardize_manifest_vars(Tall, 'STEP1 read manifest');
        fprintf('[READ] %s\n', manifestPath);
    catch ME
        if isfield(CFG,'AutoRebuildManifestIfInvalid') && CFG.AutoRebuildManifestIfInvalid
            fprintf('[WARN] Existing raw manifest is missing or invalid. Rebuilding STEP1 automatically.\n');
            fprintf('[WARN] Reason: %s\n', ME.message);
            Tall = build_manifest(CFG);
            Tall = standardize_manifest_vars(Tall, 'STEP1 rebuild_manifest');
            write_table_utf8bom(Tall, manifestPath);
            fprintf('[WRITE] %s\n', manifestPath);
        else
            rethrow(ME);
        end
    end
end

if CFG.RunStep2_CleanReslice
    fprintf('\n[STEP2] Clean and reslice con images to atlas space...\n');
    Tall = clean_reslice_manifest_conimages(Tall, CFG);
    Tall = standardize_manifest_vars(Tall, 'STEP2 clean_reslice_manifest_conimages');
    write_table_utf8bom(Tall, manifestAtlasPath);
    fprintf('[WRITE] %s\n', manifestAtlasPath);
else
    fprintf('\n[STEP2-SKIP] Read existing atlas-space manifest...\n');
    try
        assert(exist(manifestAtlasPath,'file') == 2, 'Missing atlas-space manifest CSV: %s', manifestAtlasPath);
        TallAtlas = safe_readtable(manifestAtlasPath);
        TallAtlas = standardize_manifest_vars(TallAtlas, 'STEP2 read atlas-space manifest');
        Tall = TallAtlas;
        fprintf('[READ] %s\n', manifestAtlasPath);
    catch ME
        if isfield(CFG,'AutoRebuildManifestIfInvalid') && CFG.AutoRebuildManifestIfInvalid
            fprintf('[WARN] Preparing atlas-space contrast images from the input manifest.\n');
            fprintf('[WARN] Reason: %s\n', ME.message);
            Tall = standardize_manifest_vars(Tall, 'STEP2 fallback raw manifest');
            Tall = clean_reslice_manifest_conimages(Tall, CFG);
            Tall = standardize_manifest_vars(Tall, 'STEP2 fallback clean_reslice');
            write_table_utf8bom(Tall, manifestAtlasPath);
            fprintf('[WRITE] %s\n', manifestAtlasPath);
        else
            rethrow(ME);
        end
    end
end

Tall = standardize_manifest_vars(Tall, 'Before STEP3A');
oneSampleContrastInfo = table();

if isfield(CFG,'RunStep3A_TaskOneSample') && CFG.RunStep3A_TaskOneSample
    fprintf('\n[STEP3A] Pooled one-sample task activation models...\n');
    for ci = 1:numel(CFG.ContrastSets)
        cs = CFG.ContrastSets{ci};
        Tcs = Tall(strcmp(string(Tall.ContrastSet), cs.Name), :);
        if isempty(Tcs)
            fprintf('[STEP3A-SKIP] No images for ContrastSet=%s\n', cs.Name);
            continue;
        end
        Tinfo1 = run_task_onesample_models(Tcs, cs, CFG);
        if ~isempty(Tinfo1)
            oneSampleContrastInfo = [oneSampleContrastInfo; Tinfo1]; %#ok<AGROW>
        end
    end
    oneInfoPath = fullfile(CFG.OutRoot, 'Manifest', 'SecondLevel_TaskOneSample_ContrastInfo.csv');
    if ~isempty(oneSampleContrastInfo)
        write_table_utf8bom(oneSampleContrastInfo, oneInfoPath);
        fprintf('[WRITE] %s\n', oneInfoPath);
    end
else
    fprintf('\n[STEP3A-SKIP] One-sample task models skipped.\n');
end

Tall = standardize_manifest_vars(Tall, 'Before STEP3');
allContrastInfo = table();
if CFG.RunStep3_EstimateSecondLevel || CFG.RunStep3_AddContrasts
    for ci = 1:numel(CFG.ContrastSets)
        cs = CFG.ContrastSets{ci};
        fprintf('\n============================================================\n');
        fprintf('[STEP3] Full Group x Task FF model: %s\n', cs.Name);
        fprintf('============================================================\n');

        Tcs = Tall(strcmp(string(Tall.ContrastSet), cs.Name), :);
        if isempty(Tcs)
            fprintf('[SKIP] No images for ContrastSet=%s\n', cs.Name);
            continue;
        end

        if CFG.RunStep3_EstimateSecondLevel
            [SPM, rowInfo, outDir] = run_fullmodel_ff(Tcs, cs, CFG);
        else
            outDir = fullfile(CFG.OutRoot, 'SecondLevel_FullModel', cs.Name);
            spmmat = fullfile(outDir, 'SPM.mat');
            rowInfoPath = fullfile(outDir, 'FullModel_RowInfo.csv');
            assert(exist(spmmat,'file') == 2, 'Missing SPM.mat for STEP3-SKIP: %s', spmmat);
            assert(exist(rowInfoPath,'file') == 2, 'Missing RowInfo CSV for STEP3-SKIP: %s', rowInfoPath);
            S = load(spmmat); SPM = S.SPM;
            rowInfo = readtable(rowInfoPath, 'TextType','string');
            fprintf('[STEP3-ESTIMATE-SKIP] Loaded existing SPM.mat and RowInfo for %s\n', cs.Name);
        end

        if CFG.RunStep3_AddContrasts
            Tinfo = add_fullmodel_contrasts(SPM, rowInfo, outDir, cs, CFG);
            allContrastInfo = [allContrastInfo; Tinfo]; %#ok<AGROW>
        else
            fprintf('[STEP3-CONTRAST-SKIP] Contrast creation skipped for %s\n', cs.Name);
        end
    end

    contrastInfoPath = fullfile(CFG.OutRoot, 'Manifest', 'SecondLevel_FullModel_ContrastInfo.csv');
    if ~isempty(allContrastInfo)
        write_table_utf8bom(allContrastInfo, contrastInfoPath);
        fprintf('[WRITE] %s\n', contrastInfoPath);
    end
else
    fprintf('\n[STEP3-SKIP] Skip second-level model and contrast creation.\n');
end

if CFG.RunStep4_FDROutputs
    fprintf('\n[STEP4] Voxel-FDR q < %.3f | OneSample=%s | TaskDiff=%s.\n', CFG.FDR_q, CFG.FDR_OneSample_CorrectionScope, CFG.FDR_CorrectionScope);
    run_voxel_fdr_outputs(CFG, LUT);
else
    fprintf('\n[STEP4-SKIP] Voxel-FDR output skipped.\n');
end

fprintf('\n[ALL DONE] Group activation analysis finished.\n');

end

function ContrastSets = filter_contrast_sets(ContrastSets, namesToRun)
if isempty(namesToRun)
    return;
end
namesToRun = string(namesToRun);
keep = false(size(ContrastSets));
allNames = strings(size(ContrastSets));
for i = 1:numel(ContrastSets)
    allNames(i) = string(ContrastSets{i}.Name);
    keep(i) = any(strcmp(allNames(i), namesToRun));
end
missing = setdiff(namesToRun, allNames);
if ~isempty(missing)
    error('Unknown ContrastSetsToRun: %s', strjoin(missing, ', '));
end
ContrastSets = ContrastSets(keep);
fprintf('[CONFIG] Running selected ContrastSets: %s\n', strjoin(allNames(keep), ', '));
end

function Tall = build_manifest(CFG)
rows = {};
for gi = 1:numel(CFG.Groups)
    gName = CFG.Groups{gi};
    for ri = 1:numel(CFG.Runs)
        rName = CFG.Runs{ri};
        taskName = task_from_group_run(gName, rName);

        runDir = fullfile(CFG.FirstLevelRoot, gName, rName);
        if exist(runDir,'dir') ~= 7
            fprintf('[WARN] Missing run dir: %s\n', runDir);
            continue;
        end

        subDirs = dir(fullfile(runDir, 'UP_*'));
        subDirs = subDirs([subDirs.isdir]);

        for si = 1:numel(subDirs)
            subFolder = subDirs(si).name;
            sub = erase(subFolder, 'UP_');
            subDir = fullfile(runDir, subFolder);
            spmmat = fullfile(subDir, 'SPM.mat');
            if exist(spmmat,'file') ~= 2
                fprintf('[WARN] Missing SPM.mat: %s\n', spmmat);
                continue;
            end

            try
                S = load(spmmat);
                SPM = S.SPM;
            catch ME
                fprintf('[WARN] Cannot load SPM.mat: %s\n%s\n', spmmat, ME.message);
                continue;
            end

            for ci = 1:numel(CFG.ContrastSets)
                cs = CFG.ContrastSets{ci};
                conIdx = find_contrast_index(SPM, cs.FirstLevelConName);
                if isempty(conIdx)
                    fprintf('[WARN] Contrast not found: %s | %s\n', cs.FirstLevelConName, spmmat);
                    continue;
                end

                conPath = get_contrast_path(SPM, conIdx, subDir);
                if exist(conPath,'file') ~= 2
                    fprintf('[WARN] Missing con image: %s\n', conPath);
                    continue;
                end

                rows(end+1,:) = { ...
                    char(cs.Name), char(cs.FirstLevelConName), ...
                    char(gName), gi, char(rName), ri, char(taskName), task_index(taskName, CFG), ...
                    char(sub), char(sprintf('%s_%s', gName, sub)), ...
                    conIdx, char(spmmat), char(conPath), ''}; %#ok<AGROW>
            end
        end
    end
end

Tall = cell2table(rows, 'VariableNames', { ...
    'ContrastSet','FirstLevelConName', ...
    'Group','GroupLevel','Run','RunLevel','Task','TaskLevel', ...
    'Subject','SubjectKey','FirstLevelConIndex','SPMmat','ConImage','ConImage_AtlasSpace' ...
    });

if ~isempty(Tall)
    Tall = sortrows(Tall, {'ContrastSet','GroupLevel','SubjectKey','TaskLevel'});
end

end

function T = safe_readtable(csvPath)
try
    T = readtable(csvPath, 'TextType','string', 'VariableNamingRule','preserve');
catch
    T = readtable(csvPath, 'TextType','string');
end
end

function T = standardize_manifest_vars(T, context)

required = {'ContrastSet','FirstLevelConName','Group','GroupLevel','Run','RunLevel', ...
    'Task','TaskLevel','Subject','SubjectKey','FirstLevelConIndex','SPMmat','ConImage','ConImage_AtlasSpace'};

if width(T) >= numel(required)
    vars0 = string(T.Properties.VariableNames);
    if ~any(strcmpi(strip_bom(vars0), "ContrastSet")) && height(T) >= 1
        try
            firstRow = string(table2cell(T(1,1:numel(required))));
            firstRowClean = strip_bom(firstRow);
            if any(strcmpi(firstRowClean, "ContrastSet")) && any(strcmpi(firstRowClean, "ConImage"))
                T.Properties.VariableNames(1:numel(required)) = cellstr(matlab.lang.makeValidName(firstRowClean));
                T(1,:) = [];
            end
        catch
        end
    end
end

vars = string(T.Properties.VariableNames);
varsClean = strip_bom(strtrim(vars));
normVars = normalize_header(varsClean);

for k = 1:numel(required)
    req = string(required{k});
    reqNorm = normalize_header(req);

    idx = find(normVars == reqNorm, 1, 'first');

    if isempty(idx)
        idx = find(endsWith(normVars, reqNorm), 1, 'first');
    end

    if ~isempty(idx)
        T.Properties.VariableNames{idx} = char(req);
    end
end

missing = setdiff(required, T.Properties.VariableNames, 'stable');
if ~isempty(missing)
    fprintf('\n[ERROR] Manifest variables are incomplete at %s.\n', context);
    fprintf('[ERROR] Existing variables are:\n');
    disp(T.Properties.VariableNames(:));
    fprintf('[ERROR] Missing variables are:\n');
    disp(missing(:));
    error('Manifest table does not contain the required variables.');
end
end

function s = strip_bom(s)
s = string(s);
s = replace(s, char(65279), "");
s = regexprep(s, '^\xEF\xBB\xBF', '');
s = strtrim(s);
end

function s = normalize_header(s)
s = strip_bom(string(s));
s = lower(regexprep(s, '[^A-Za-z0-9]', ''));
end

function conIdx = find_contrast_index(SPM, conName)
conIdx = [];
if ~isfield(SPM,'xCon') || isempty(SPM.xCon), return; end
names = string({SPM.xCon.name});
hit = find(strcmp(names, string(conName)), 1, 'first');
if isempty(hit)
    names2 = erase(names, " ");
    target2 = erase(string(conName), " ");
    hit = find(strcmp(names2, target2), 1, 'first');
end
if ~isempty(hit), conIdx = hit; end
end

function conPath = get_contrast_path(SPM, conIdx, subDir)
conPath = '';
try
    if isfield(SPM.xCon(conIdx), 'Vcon') && ~isempty(SPM.xCon(conIdx).Vcon)
        fname = SPM.xCon(conIdx).Vcon.fname;
        if iscell(fname), fname = fname{1}; end
        if ischar(fname) || isstring(fname)
            fname = char(fname);
            if is_absolute_path(fname)
                conPath = fname;
            else
                conPath = fullfile(subDir, fname);
            end
        end
    end
catch
end

if isempty(conPath)
    conPath = fullfile(subDir, sprintf('con_%04d.nii', conIdx));
end
end

function tf = is_absolute_path(p)
p = char(p);
tf = numel(p) >= 2 && p(2) == ':';
if ~tf
    tf = startsWith(p, filesep);
end
end

function Tall = clean_reslice_manifest_conimages(Tall, CFG)
outRoot = fullfile(CFG.OutRoot, 'ConAtlasSpace');

for i = 1:height(Tall)
    src = char(Tall.ConImage{i});
    cs = char(Tall.ContrastSet{i});
    gName = char(Tall.Group{i});
    rName = char(Tall.Run{i});
    taskName = char(Tall.Task{i});
    sub = char(Tall.Subject{i});

    outDir = fullfile(outRoot, cs, gName, rName);
    mkdir_if_needed(outDir);

    outName = sprintf('%s_%s_%s_%s_%s_atlas.nii', cs, gName, rName, taskName, sub);
    outName = regexprep(outName, '[^\w\-\.\(\)]', '_');
    outPath = fullfile(outDir, outName);

    if CFG.SkipExistingAtlasConImages && exist(outPath,'file') == 2
        current = spm_vol(outPath); target = spm_vol(CFG.AtlasNii);
        if isequal(current.dim,target.dim) && max(abs(current.mat(:)-target.mat(:)))<1e-4
            Tall.ConImage_AtlasSpace{i} = outPath;
            continue;
        end
    end

    if CFG.CleanAndResliceConImages
        outPath = clean_and_reslice_one(src, outPath, CFG);
    else
        outPath = src;
    end

    Tall.ConImage_AtlasSpace{i} = outPath;
end
end

function outPath = clean_and_reslice_one(src, outPath, CFG)
mkdir_if_needed(fileparts(outPath));

Vatlas = spm_vol(CFG.AtlasNii);
Vsrc = spm_vol(src);
Y = spm_read_vols(Vsrc);
Y(~isfinite(Y)) = 0;

tmpDir = fullfile(fileparts(outPath), '_clean_tmp');
mkdir_if_needed(tmpDir);

[~, base, ~] = fileparts(src);
cleanTmp = fullfile(tmpDir, ['clean_' base '.nii']);

Vclean = Vsrc;
Vclean.fname = cleanTmp;
Vclean.dt = [16 0];
spm_write_vol(Vclean, Y);

sameGrid = isequal(Vclean.dim, Vatlas.dim) && max(abs(Vclean.mat(:) - Vatlas.mat(:))) < 1e-4;

if sameGrid
    copyfile(cleanTmp, outPath);
else
    flags = struct('interp', CFG.ResliceInterp, 'wrap', [0 0 0], 'mask', 0, 'mean', 0, 'which', 1, 'prefix', 'r');
    P = char(CFG.AtlasNii, cleanTmp);
    spm_reslice(P, flags);

    rfile = fullfile(tmpDir, ['rclean_' base '.nii']);
    if exist(rfile,'file') ~= 2
        error('Resliced file not found: %s', rfile);
    end
    movefile(rfile, outPath, 'f');
end

Vout = spm_vol(outPath);
if ~isequal(Vout.dim, Vatlas.dim) || max(abs(Vout.mat(:) - Vatlas.mat(:))) > 1e-4
    error('Output con image is not in atlas grid: %s', outPath);
end
end

function TinfoAll = run_task_onesample_models(Tcs, cs, CFG)
% A task-level contrast image contributes once per participant to an
% intercept-only model. Positive activation tests its pooled mean against zero.
if CFG.RequireCompleteSubject
    keys = unique(string(Tcs.SubjectKey),'stable'); included = false(height(Tcs),1);
    for si = 1:numel(keys)
        rows = string(Tcs.SubjectKey)==keys(si);
        if all(ismember(string(CFG.Tasks),unique(string(Tcs.Task(rows)))))
            included = included | rows;
        end
    end
    Tcs = Tcs(included,:);
end
TinfoAll = table();
for ti = 1:numel(CFG.Tasks)
    taskName = CFG.Tasks{ti};
    Ttask = Tcs(string(Tcs.Task)==taskName,:);
    assert(height(Ttask)>=2,'At least two participant contrasts are required for %s.',taskName);
    assert(numel(unique(string(Ttask.SubjectKey)))==height(Ttask), ...
        'Each participant must contribute exactly one contrast to %s.',taskName);
    scans = cellstr(string(Ttask.ConImage_AtlasSpace));
    for k = 1:numel(scans)
        file = regexprep(scans{k},',\d+$','');
        assert(isfile(file),'Missing subject contrast: %s',file);
        if isempty(regexp(scans{k},',\d+$','once')), scans{k} = [scans{k} ',1']; end
    end
    outDir = fullfile(CFG.OutRoot,'SecondLevel_TaskOneSample',cs.Name,taskName);
    spmmat = fullfile(outDir,'SPM.mat');
    if isfolder(outDir) && CFG.OverwriteTaskOneSample, rmdir(outDir,'s'); end
    mkdir_if_needed(outDir);
    if ~isfile(spmmat)
        clear matlabbatch
        matlabbatch{1}.spm.stats.factorial_design.dir = {outDir};
        matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = scans(:);
        matlabbatch{1}.spm.stats.factorial_design.cov = struct('c',{},'cname',{},'iCFI',{},'iCC',{});
        matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
        if CFG.UseExplicitROIMaskInSPMModel
            matlabbatch{1}.spm.stats.factorial_design.masking.em = {[CFG.UnionMaskNii ',1']};
        else
            matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
        end
        matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep( ...
            'Factorial design specification: SPM.mat File', ...
            substruct('.','val','{}',{1},'.','val','{}',{1},'.','val','{}',{1}), ...
            substruct('.','spmmat'));
        matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
        save(fullfile(outDir,'SPM_batch_one_sample.mat'),'matlabbatch');
        spm_jobman('run',matlabbatch);
    end
    loaded = load(spmmat,'SPM'); SPM = loaded.SPM;
    assert(size(SPM.xX.X,2)==1,'Task activation requires an intercept-only one-sample design.');
    if CFG.RecreateContrasts && isfield(SPM,'xCon'), SPM = rmfield(SPM,'xCon'); end
    if ~isfield(SPM,'xCon'), SPM.xCon = struct([]); end
    cname = sprintf('%s_OneSample_%s_gt_0',cs.Name,taskName);
    conIndex = find_contrast_index(SPM,cname);
    if isempty(conIndex)
        SPM.xCon = [SPM.xCon spm_FcUtil('Set',cname,'T','c',1,SPM.xX.xKXs)];
        conIndex = numel(SPM.xCon);
    end
    save(spmmat,'SPM');
    SPM = spm_contrasts(SPM,conIndex);
    save(spmmat,'SPM');
    Tinfo = cell2table({cs.Name,cname,'T','OneSample',taskName,'',outDir,conIndex}, ...
        'VariableNames',{'ContrastSet','ContrastName','StatType','AnalysisType', ...
        'Task','Group','ModelDir','LocalContrastOrder'});
    write_table_utf8bom(Tinfo,fullfile(outDir,'OneSample_ContrastInfo.csv'));
    write_table_utf8bom(Ttask,fullfile(outDir,'OneSample_InputManifest.csv'));
    TinfoAll = [TinfoAll; Tinfo]; %#ok<AGROW>
end
end

function [SPM, rowInfo, outDir] = run_fullmodel_ff(Tcs, cs, CFG)
% Subject intercepts model repeated observations; GroupTaskCell represents
% all 16 Latin-square group-by-task cells. Marginal task differences are
% estimable contrasts, with equal weight assigned to each group.
outDir = fullfile(CFG.OutRoot, 'SecondLevel_FullModel', cs.Name);

if exist(outDir,'dir') == 7 && isfield(CFG,'OverwriteSecondLevel') && CFG.OverwriteSecondLevel
    fprintf('[STEP3-OVERWRITE] Removing existing full-model directory: %s\n', outDir);
    try
        rmdir(outDir, 's');
    catch ME
        error('Cannot remove existing full-model output directory: %s\n%s', outDir, ME.message);
    end
end
mkdir_if_needed(outDir);

spmmat = fullfile(outDir, 'SPM.mat');

subjects = unique(string(Tcs.SubjectKey), 'stable');
keep = false(numel(subjects),1);

for si = 1:numel(subjects)
    Ts = Tcs(strcmp(string(Tcs.SubjectKey), subjects(si)), :);
    hasAll = true;
    for ti = 1:numel(CFG.Tasks)
        hasAll = hasAll && any(strcmp(string(Ts.Task), CFG.Tasks{ti}));
    end
    keep(si) = hasAll || ~CFG.RequireCompleteSubject;
end

subjects = subjects(keep);
if isempty(subjects)
    error('No subjects available for full model: %s', cs.Name);
end

fprintf('[INFO] %s | complete-case subjects = %d\n', cs.Name, numel(subjects));

rowRows = {};
scanRows = {};
subLevels = [];
groupLevels = [];
taskLevels = [];
cellLevels = [];
cellNames = {};

for si = 1:numel(subjects)
    skey = subjects(si);
    Ts = Tcs(strcmp(string(Tcs.SubjectKey), skey), :);

    glev = unique(double(Ts.GroupLevel));
    if numel(glev) ~= 1
        error('SubjectKey has multiple groups: %s', skey);
    end

    for ti = 1:numel(CFG.Tasks)
        taskName = CFG.Tasks{ti};
        idx = find(strcmp(string(Ts.Task), taskName), 1, 'first');
        if isempty(idx)
            if CFG.RequireCompleteSubject
                error('Missing task %s for subject %s', taskName, skey);
            else
                continue;
            end
        end

        img = char(Ts.ConImage_AtlasSpace{idx});
        cellLevel = (glev(1)-1) * numel(CFG.Tasks) + ti;
        cellName = sprintf('G%d_%s', glev(1), taskName);

        scanRows{end+1,1} = [img ',1']; %#ok<AGROW>
        subLevels(end+1,1) = si; %#ok<AGROW>
        groupLevels(end+1,1) = glev(1); %#ok<AGROW>
        taskLevels(end+1,1) = ti; %#ok<AGROW>
        cellLevels(end+1,1) = cellLevel; %#ok<AGROW>
        cellNames{end+1,1} = cellName; %#ok<AGROW>

        rowRows(end+1,:) = { ...
            numel(scanRows), char(cs.Name), char(skey), char(Ts.Subject{idx}), ...
            char(Ts.Group{idx}), glev(1), char(taskName), ti, ...
            cellLevel, char(cellName), ...
            char(Ts.Run{idx}), double(Ts.RunLevel(idx)), img}; %#ok<AGROW>
    end
end

rowInfo = cell2table(rowRows, 'VariableNames', { ...
    'ScanIndex','ContrastSet','SubjectKey','Subject','Group','GroupLevel','Task','TaskLevel', ...
    'CellLevel','CellName','Run','RunLevel','ConImage_AtlasSpace' ...
    });

rowInfoPath = fullfile(outDir, 'FullModel_RowInfo.csv');
write_table_utf8bom(rowInfo, rowInfoPath);
fprintf('[WRITE] %s\n', rowInfoPath);

if isfield(CFG, 'RunPreflightVoxelQC') && CFG.RunPreflightVoxelQC
    preflight_model_inputs(scanRows, CFG, cs.Name, outDir);
end

if exist(spmmat,'file') == 2 && ~CFG.OverwriteSecondLevel
    S = load(spmmat);
    SPM = S.SPM;
    return;
end

clear matlabbatch;
matlabbatch{1}.spm.stats.factorial_design.dir = {outDir};

matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).name = 'Subject';
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).variance = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).ancova = 0;

matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).name = 'GroupTaskCell';
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).dept = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).variance = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).ancova = 0;

matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{1}.fmain.fnum = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{2}.fmain.fnum = 2;

for si = 1:numel(subjects)
    rows = find(subLevels == si);
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(si).scans = scanRows(rows);
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(si).conds = [subLevels(rows), cellLevels(rows)];
end

matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;

if isfield(CFG, 'UseExplicitROIMaskInSPMModel') && CFG.UseExplicitROIMaskInSPMModel
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {[CFG.UnionMaskNii ',1']};
    fprintf('[MASK] SPM estimation uses explicit ROI mask: %s\n', CFG.UnionMaskNii);
else
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
    fprintf('[MASK] SPM estimation uses no explicit ROI mask; FDR outputs remain ROI-restricted.\n');
end

matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep( ...
    'Factorial design specification: SPM.mat File', ...
    substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

batchPath = fullfile(outDir, sprintf('SPM_batch_design_estimate_%s.mat', safe_filename(cs.Name)));
save(batchPath, 'matlabbatch');
fprintf('[WRITE] %s\n', batchPath);

try
    spm_jobman('run', matlabbatch);
catch ME
    fprintf(2, '\n[SPM-ERROR] Step3 design/estimation failed for %s\n', cs.Name);
    fprintf(2, '[SPM-ERROR] Batch saved at: %s\n', batchPath);
    fprintf(2, '[SPM-ERROR] Message: %s\n', ME.message);
    if exist(fullfile(outDir,'SPM.mat'), 'file') == 2
        try
            Sfail = load(fullfile(outDir,'SPM.mat'));
            if isfield(Sfail,'SPM') && isfield(Sfail.SPM,'xX') && isfield(Sfail.SPM.xX,'X')
                Xfail = Sfail.SPM.xX.X;
                fprintf(2, '[SPM-ERROR] Design rows=%d cols=%d rank=%d\n', size(Xfail,1), size(Xfail,2), rank(full(Xfail)));
            end
        catch
        end
    end
    rethrow(ME);
end

S = load(spmmat);
SPM = S.SPM;
end

function Tinfo = add_fullmodel_contrasts(SPM, rowInfo, outDir, cs, CFG)
spmmat = fullfile(outDir, 'SPM.mat');
if nargin < 1 || isempty(SPM)
    S = load(spmmat);
    SPM = S.SPM;
end

nG = numel(CFG.Groups);
nT = numel(CFG.Tasks);
nCell = nG * nT;
ncol = size(SPM.xX.X,2);

cellColsVec = find_tail_hypothesis_cols_ff(SPM, nCell, outDir, 'GroupTaskCell');
cellCols = reshape(cellColsVec, [nT, nG])';  % rows=Group, columns=Task if cells were entered G-major/T-minor

if isfield(CFG,'DebugContrastColumns') && CFG.DebugContrastColumns
    fp = fullfile(outDir, sprintf('ContrastColumns_%s.txt', safe_filename(cs.Name)));
    fid = fopen(fp,'w');
    fprintf(fid, 'Cell columns in G-major/T-minor order:\n');
    for gi = 1:nG
        for ti = 1:nT
            fprintf(fid, 'G%d_%s -> col %d : %s\n', gi, CFG.Tasks{ti}, cellCols(gi,ti), char(string(SPM.xX.name{cellCols(gi,ti)})));
        end
    end
    fclose(fid);
end

if isfield(SPM,'xCon') && CFG.RecreateContrasts
    SPM = rmfield(SPM, 'xCon');
end
if ~isfield(SPM,'xCon') || isempty(SPM.xCon)
    SPM.xCon = struct([]);
end

infoRows = {};
conCounter = 0;

if CFG.EqualGroupWeight
    gWeights = ones(nG,1) ./ nG;
else
    gWeights = zeros(nG,1);
    for gi = 1:nG
        gWeights(gi) = sum(rowInfo.GroupLevel == gi);
    end
    gWeights = gWeights ./ sum(gWeights);
end

taskRows = zeros(nT, ncol);
for ti = 1:nT
    w = zeros(1,ncol);
    for gi = 1:nG
        w(cellCols(gi,ti)) = gWeights(gi);
    end
    taskRows(ti,:) = w;
end

Wf = zeros(3, ncol);
Wf(1,:) = taskRows(1,:) - taskRows(2,:); % Safety - Vitality
Wf(2,:) = taskRows(1,:) - taskRows(3,:); % Safety - Beauty
Wf(3,:) = taskRows(1,:) - taskRows(4,:); % Safety - Imageability
cname = sprintf('%s_TaskMainEffect_F', cs.Name);
SPM.xCon = [SPM.xCon spm_FcUtil('Set', cname, 'F', 'c', Wf', SPM.xX.xKXs)];
conCounter = conCounter + 1;
infoRows(end+1,:) = {cs.Name, cname, 'F', 'TaskMainEffect', '', '', '', conCounter}; %#ok<AGROW>

for pi = 1:size(CFG.PairList,1)
    A = CFG.PairList{pi,1};
    B = CFG.PairList{pi,2};
    ia = task_index(A, CFG);
    ib = task_index(B, CFG);

    w = taskRows(ia,:) - taskRows(ib,:);
    cname = sprintf('%s_%s_gt_%s', cs.Name, A, B);
    SPM.xCon = [SPM.xCon spm_FcUtil('Set', cname, 'T', 'c', w', SPM.xX.xKXs)];
    conCounter = conCounter + 1;
    infoRows(end+1,:) = {cs.Name, cname, 'T', 'Pairwise', '', A, B, conCounter}; %#ok<AGROW>

    if ~isfield(CFG, 'WriteBidirectionalPairwise') || CFG.WriteBidirectionalPairwise
        w = taskRows(ib,:) - taskRows(ia,:);
        cname = sprintf('%s_%s_gt_%s', cs.Name, B, A);
        SPM.xCon = [SPM.xCon spm_FcUtil('Set', cname, 'T', 'c', w', SPM.xX.xKXs)];
        conCounter = conCounter + 1;
        infoRows(end+1,:) = {cs.Name, cname, 'T', 'Pairwise', '', B, A, conCounter}; %#ok<AGROW>
    end
end

if isfield(CFG, 'AddPairwiseEqualityF') && CFG.AddPairwiseEqualityF
    for pi = 1:size(CFG.PairList,1)
        A = CFG.PairList{pi,1};
        B = CFG.PairList{pi,2};
        ia = task_index(A, CFG);
        ib = task_index(B, CFG);

        w = taskRows(ia,:) - taskRows(ib,:);
        cname = sprintf('%s_%s_vs_%s_%s', cs.Name, A, B, CFG.PairwiseEqualityF_Suffix);
        SPM.xCon = [SPM.xCon spm_FcUtil('Set', cname, 'F', 'c', w', SPM.xX.xKXs)];
        conCounter = conCounter + 1;
        infoRows(end+1,:) = {cs.Name, cname, 'F', 'PairwiseEqualityF', '', A, B, conCounter}; %#ok<AGROW>
    end
end

if CFG.AddInteractionF
    refG = 1;
    refT = 1;
    Wint = zeros((nG-1)*(nT-1), ncol);
    rr = 0;
    for gi = 2:nG
        for ti = 2:nT
            rr = rr + 1;
            Wint(rr, cellCols(gi,ti))   =  1;
            Wint(rr, cellCols(gi,refT)) = -1;
            Wint(rr, cellCols(refG,ti)) = -1;
            Wint(rr, cellCols(refG,refT)) = 1;
        end
    end
    cname = sprintf('%s_TaskXGroup_Interaction_F', cs.Name);
    SPM.xCon = [SPM.xCon spm_FcUtil('Set', cname, 'F', 'c', Wint', SPM.xX.xKXs)];
    conCounter = conCounter + 1;
    infoRows(end+1,:) = {cs.Name, cname, 'F', 'TaskXGroupInteraction', '', '', '', conCounter}; %#ok<AGROW>
end

save(spmmat, 'SPM');

SPM = spm_contrasts(SPM, 1:numel(SPM.xCon));
save(spmmat, 'SPM');

Tinfo = cell2table(infoRows, 'VariableNames', { ...
    'ContrastSet','ContrastName','StatType','AnalysisType','Task','TaskA','TaskB','LocalContrastOrder' ...
    });
write_table_utf8bom(Tinfo, fullfile(outDir, 'FullModel_ContrastInfo.csv'));
fprintf('[WRITE] %s\n', fullfile(outDir, 'FullModel_ContrastInfo.csv'));
end

function cols = find_tail_hypothesis_cols_ff(SPM, nTail, outDir, effectName)
if isfield(SPM.xX,'iH') && ~isempty(SPM.xX.iH)
    iH = SPM.xX.iH(:)';
    if numel(iH) >= nTail
        cols = iH(end-nTail+1:end);
        return;
    end
end

fp = fullfile(outDir, sprintf('FF_design_debug_%s.txt', safe_filename(effectName)));
fid = fopen(fp,'w');
if fid > 0
    xnames = SPM.xX.name(:);
    for i = 1:numel(xnames)
        fprintf(fid, '%4d : %s\n', i, char(string(xnames{i})));
    end
    fclose(fid);
end
error('Could not identify trailing FF columns for %s. See %s', effectName, fp);
end

function run_voxel_fdr_outputs(CFG, LUT)
Vroi = spm_vol(CFG.AtlasNii);
Yroi = round(spm_read_vols(Vroi));
roiMask = Yroi > 0 & isfinite(Yroi);

allVoxTables = {};
allRoiTables = {};

for ci = 1:numel(CFG.ContrastSets)
    cs = CFG.ContrastSets{ci};

    modelDirs = {};

    if isfield(CFG,'FDR_IncludeTaskOneSample') && CFG.FDR_IncludeTaskOneSample
        oneRoot = fullfile(CFG.OutRoot, 'SecondLevel_TaskOneSample', cs.Name);
        if exist(oneRoot,'dir') == 7
            for ti = 1:numel(CFG.Tasks)
                d = fullfile(oneRoot, CFG.Tasks{ti});
                if exist(fullfile(d,'SPM.mat'), 'file') == 2
                    modelDirs{end+1,1} = d; %#ok<AGROW>
                end
            end
        end
    end

    fullDir = fullfile(CFG.OutRoot, 'SecondLevel_FullModel', cs.Name);
    if exist(fullfile(fullDir,'SPM.mat'), 'file') == 2
        modelDirs{end+1,1} = fullDir; %#ok<AGROW>
    end

    if isempty(modelDirs)
        fprintf('[FDR-SKIP] No SPM.mat found for ContrastSet=%s\n', cs.Name);
        continue;
    end

    for mi = 1:numel(modelDirs)
        modelDir = modelDirs{mi};
        spmmat = fullfile(modelDir, 'SPM.mat');
        S = load(spmmat);
        SPM = S.SPM;

        if ~isfield(SPM,'xCon') || isempty(SPM.xCon)
            fprintf('[FDR-SKIP] No contrasts in: %s\n', spmmat);
            continue;
        end

        for xi = 1:numel(SPM.xCon)
            xcon = SPM.xCon(xi);
            cname = char(string(xcon.name));
            analysisType = classify_contrast_name(cname);

            if strcmpi(analysisType, 'OneSample') && ~(isfield(CFG,'FDR_IncludeTaskOneSample') && CFG.FDR_IncludeTaskOneSample)
                continue;
            end
            if strcmpi(analysisType, 'TaskXGroupInteraction') && ~(isfield(CFG,'FDR_IncludeInteractionF') && CFG.FDR_IncludeInteractionF)
                continue;
            end

            if isempty(xcon.Vspm) || ~isfield(xcon.Vspm,'fname')
                fprintf('[FDR-SKIP] No Vspm for contrast: %s\n', cname);
                continue;
            end

            statPath = fullfile(modelDir, xcon.Vspm.fname);
            if exist(statPath,'file') ~= 2
                fprintf('[FDR-SKIP] Missing stat map: %s\n', statPath);
                continue;
            end

            fdrScope = fdr_scope_for_analysis(analysisType, CFG);
            fprintf('[FDR-%s] %s | %s\n', fdrScope, cs.Name, cname);
            [Tvox, Troi] = compute_and_write_fdr_one_contrast(SPM, xi, statPath, Yroi, roiMask, LUT, CFG);

            if ~isempty(Tvox)
                allVoxTables{end+1,1} = Tvox; %#ok<AGROW>
            end
            if ~isempty(Troi)
                allRoiTables{end+1,1} = Troi; %#ok<AGROW>
            end
        end
    end
end

if ~isempty(allVoxTables)
    Tallvox = vertcat(allVoxTables{:});
    outAll = fullfile(CFG.OutRoot, 'FDR_VoxelTables', sprintf('ALL_FDR_sigvoxels_q%03d.csv', round(CFG.FDR_q*1000)));
    write_table_utf8bom(Tallvox, outAll);
    fprintf('[WRITE] %s\n', outAll);
end

if ~isempty(allRoiTables)
    Tallroi = vertcat(allRoiTables{:});
    outAll = fullfile(CFG.OutRoot, 'FDR_ROISummary', sprintf('ALL_FDR_roi_summary_q%03d.csv', round(CFG.FDR_q*1000)));
    write_table_utf8bom(Tallroi, outAll);
    fprintf('[WRITE] %s\n', outAll);
end
end

function [Tvox, Troi] = compute_and_write_fdr_one_contrast(SPM, conIdx, statPath, Yroi, roiMask, LUT, CFG)
xcon = SPM.xCon(conIdx);
cname = char(string(xcon.name));
statType = xcon.STAT;
analysisType = classify_contrast_name(cname);

Vstat = spm_vol(statPath);
Ystat = spm_read_vols(Vstat);

validAny = roiMask & isfinite(Ystat);
if ~any(validAny(:))
    Tvox = table();
    Troi = table();
    return;
end

switch upper(statType)
    case 'T'
        df1 = NaN;
        df2 = SPM.xX.erdf;
    case 'F'
        df1 = xcon.eidf;
        df2 = SPM.xX.erdf;
    otherwise
        error('Unsupported STAT type: %s', statType);
end

fdrScope = fdr_scope_for_analysis(analysisType, CFG);
if strcmpi(fdrScope, 'Union')
    [pMap, qMap, sigMask] = compute_pmap_qmap_union(Ystat, roiMask, statType, df1, df2, CFG, analysisType);
else
    [pMap, qMap, sigMask] = compute_pmap_qmap_roiwise(Ystat, Yroi, roiMask, LUT, statType, df1, df2, CFG, analysisType);
end

safeName = safe_filename(cname);
mapDir = fullfile(CFG.OutRoot, 'FDR_Maps');
mkdir_if_needed(mapDir);

if isfield(CFG, 'WriteStatMap') && CFG.WriteStatMap
    Vs = Vstat; Vs.fname = fullfile(mapDir, [safeName '_statmap.nii']); Vs.dt = [16 0];
    spm_write_vol(Vs, Ystat);
end

if isfield(CFG, 'WriteThresholdedStatMap') && CFG.WriteThresholdedStatMap
    statThrMap = zeros(size(Ystat));
    statThrMap(sigMask) = Ystat(sigMask);
    Vst = Vstat;
    Vst.fname = fullfile(mapDir, sprintf('%s_statmap_q%03d.nii', safeName, round(CFG.FDR_q*1000)));
    Vst.dt = [16 0];
    spm_write_vol(Vst, statThrMap);
end

if CFG.WritePmapQmapSigmask
    Vp = Vstat; Vp.fname = fullfile(mapDir, [safeName '_pmap.nii']); Vp.dt = [16 0];
    Vq = Vstat; Vq.fname = fullfile(mapDir, [safeName '_qmap.nii']); Vq.dt = [16 0];
    Vm = Vstat; Vm.fname = fullfile(mapDir, sprintf('%s_sigmask_q%03d.nii', safeName, round(CFG.FDR_q*1000))); Vm.dt = [2 0];
    spm_write_vol(Vp, pMap);
    spm_write_vol(Vq, qMap);
    spm_write_vol(Vm, uint8(sigMask));
end

sigIdx = find(sigMask);
if isfinite(CFG.MaxVoxelRowsPerContrast) && numel(sigIdx) > CFG.MaxVoxelRowsPerContrast
    qAll = qMap(sigIdx);
    [~, ord] = sort(qAll, 'ascend');
    sigIdx = sigIdx(ord(1:CFG.MaxVoxelRowsPerContrast));
end

Tvox = build_sigvoxel_table(sigIdx, Ystat, pMap, qMap, Yroi, Vstat, LUT, xcon, analysisType, CFG);
Troi = build_roi_summary_table(sigMask, Ystat, pMap, qMap, Yroi, Vstat, LUT, xcon, analysisType, CFG);

if CFG.WriteVoxelCSV
    voxelDir = fullfile(CFG.OutRoot, 'FDR_VoxelTables');
    mkdir_if_needed(voxelDir);
    outVox = fullfile(voxelDir, sprintf('%s_sigvoxels_q%03d.csv', safeName, round(CFG.FDR_q*1000)));
    write_table_utf8bom(Tvox, outVox);
    fprintf('[WRITE] %s | nSigVox=%d\n', outVox, height(Tvox));
end

if CFG.WriteROISummaryCSV
    roiDir = fullfile(CFG.OutRoot, 'FDR_ROISummary');
    mkdir_if_needed(roiDir);
    outRoi = fullfile(roiDir, sprintf('%s_roi_summary_q%03d.csv', safeName, round(CFG.FDR_q*1000)));
    write_table_utf8bom(Troi, outRoi);
    fprintf('[WRITE] %s\n', outRoi);
end

if isfield(CFG,'WriteFDRByROIFolders') && CFG.WriteFDRByROIFolders && should_write_fdr_byroi_files(analysisType)
    write_fdr_byroi_files(Tvox, Troi, safeName, cname, analysisType, LUT, CFG);
end
end

function [pMap, qMap, sigMask] = compute_pmap_qmap_union(Ystat, roiMask, statType, df1, df2, CFG, analysisType)
if nargin < 7
    analysisType = '';
end
valid = roiMask & isfinite(Ystat);
statVec = Ystat(valid);
pVec = p_from_stat(statVec, statType, df1, df2, CFG, analysisType);
qVec = bh_fdr_q(pVec);
pMap = nan(size(Ystat));
qMap = nan(size(Ystat));
pMap(valid) = pVec;
qMap(valid) = qVec;
sigMask = false(size(Ystat));
sigMask(valid) = qVec <= CFG.FDR_q;
end

% Benjamini-Hochberg correction is evaluated independently within each ROI.
function [pMap, qMap, sigMask] = compute_pmap_qmap_roiwise(Ystat, Yroi, roiMask, LUT, statType, df1, df2, CFG, analysisType)
if nargin < 9
    analysisType = '';
end
pMap = nan(size(Ystat));
qMap = nan(size(Ystat));
sigMask = false(size(Ystat));

for r = 1:numel(LUT.ROI_ID)
    id = LUT.ROI_ID(r);
    valid = (Yroi == id) & roiMask & isfinite(Ystat);
    if ~any(valid(:))
        continue;
    end
    statVec = Ystat(valid);
    pVec = p_from_stat(statVec, statType, df1, df2, CFG, analysisType);
    qVec = bh_fdr_q(pVec);
    pMap(valid) = pVec;
    qMap(valid) = qVec;
    sigMask(valid) = qVec <= CFG.FDR_q;
end
end

function pVec = p_from_stat(statVec, statType, df1, df2, CFG, analysisType)
switch upper(statType)
    case 'T'
        if strcmpi(analysisType, 'Pairwise') && ...
                isfield(CFG,'PairwiseT_TwoSidedFDR') && CFG.PairwiseT_TwoSidedFDR
            pVec = 2 .* (1 - spm_Tcdf(abs(statVec), df2));
        else
            pVec = 1 - spm_Tcdf(statVec, df2);
        end
    case 'F'
        pVec = 1 - spm_Fcdf(statVec, df1, df2);
    otherwise
        error('Unsupported STAT type: %s', statType);
end
pVec(~isfinite(pVec)) = 1;
pVec = max(min(pVec, 1), 0);
end

function tf = should_write_fdr_byroi_files(analysisType)
tf = any(strcmpi(analysisType, {'TaskMainEffect','Pairwise','PairwiseEqualityF'}));
end

function fdrScope = fdr_scope_for_analysis(analysisType, CFG)
if strcmpi(analysisType, 'OneSample')
    if isfield(CFG,'FDR_OneSample_CorrectionScope') && ~isempty(CFG.FDR_OneSample_CorrectionScope)
        fdrScope = CFG.FDR_OneSample_CorrectionScope;
    else
        fdrScope = 'Union';
    end
else
    if isfield(CFG,'FDR_CorrectionScope') && ~isempty(CFG.FDR_CorrectionScope)
        fdrScope = CFG.FDR_CorrectionScope;
    else
        fdrScope = 'ROI';
    end
end
end

function write_fdr_byroi_files(Tvox, Troi, safeName, contrastName, analysisType, LUT, CFG)
[analysisFolder, contrastFolder] = contrast_subdir_parts(contrastName, analysisType, CFG);

for r = 1:numel(LUT.ROI_ID)
    id = LUT.ROI_ID(r);
    roiFolder = sprintf('%02d_%s', id, safe_filename(LUT.ROI_Abbrev{r}));

    if CFG.WriteVoxelCSV
        voxelDir = fullfile(CFG.OutRoot, 'FDR_VoxelTables', 'ByROI', roiFolder, analysisFolder, contrastFolder);
        mkdir_if_needed(voxelDir);
        if ~isempty(Tvox) && height(Tvox) > 0
            Tv = Tvox(Tvox.Merged_ROI_ID == id, :);
        else
            Tv = empty_voxel_table();
        end
        outVox = fullfile(voxelDir, sprintf('%s_%s_sigvoxels_q%03d.csv', safeName, roiFolder, round(CFG.FDR_q*1000)));
        write_table_utf8bom(Tv, outVox);
    end

    if CFG.WriteROISummaryCSV
        roiDir = fullfile(CFG.OutRoot, 'FDR_ROISummary', 'ByROI', roiFolder, analysisFolder, contrastFolder);
        mkdir_if_needed(roiDir);
        if ~isempty(Troi) && height(Troi) > 0
            Tr = Troi(Troi.Merged_ROI_ID == id, :);
        else
            Tr = empty_roi_summary_table(LUT, r, contrastName, analysisType);
        end
        outRoi = fullfile(roiDir, sprintf('%s_%s_roi_summary_q%03d.csv', safeName, roiFolder, round(CFG.FDR_q*1000)));
        write_table_utf8bom(Tr, outRoi);
    end
end
end

function [analysisFolder, contrastFolder] = contrast_subdir_parts(contrastName, analysisType, CFG)
contrastName = char(string(contrastName));
switch lower(analysisType)
    case 'taskmaineffect'
        analysisFolder = 'TaskMainEffect';
        contrastFolder = 'TaskMainEffect_F';
    case 'pairwise'
        analysisFolder = 'Pairwise';
        [A, B] = parse_pairwise_gt_name(contrastName, CFG);
        contrastFolder = sprintf('%s_gt_%s', A, B);
    case 'pairwiseequalityf'
        analysisFolder = 'PairwiseEqualityF';
        [A, B] = parse_pairwise_equality_name(contrastName, CFG);
        contrastFolder = sprintf('%s_vs_%s', A, B);
    otherwise
        analysisFolder = safe_filename(analysisType);
        contrastFolder = safe_filename(contrastName);
end
analysisFolder = safe_filename(analysisFolder);
contrastFolder = safe_filename(contrastFolder);
end

function [taskA, taskB] = parse_pairwise_gt_name(cname, CFG)
taskA = 'TaskA';
taskB = 'TaskB';
cname = char(string(cname));
for pi = 1:size(CFG.PairList,1)
    A = CFG.PairList{pi,1};
    B = CFG.PairList{pi,2};
    tokenAB = sprintf('_%s_gt_%s', A, B);
    tokenBA = sprintf('_%s_gt_%s', B, A);
    if contains(cname, tokenAB)
        taskA = A; taskB = B; return;
    elseif contains(cname, tokenBA)
        taskA = B; taskB = A; return;
    end
end
parts = split(string(cname), '_');
idx = find(parts == "gt", 1, 'last');
if ~isempty(idx) && idx > 1 && idx < numel(parts)
    taskA = char(parts(idx-1));
    taskB = char(parts(idx+1));
end
end

function T = empty_roi_summary_table(LUT, lutRow, contrastName, analysisType)
ContrastSet = {parse_contrastset(contrastName)};
ContrastName = {char(string(contrastName))};
AnalysisType = {char(string(analysisType))};
StatType = {''};
Merged_ROI_ID = LUT.ROI_ID(lutRow);
ROI_Abbrev = LUT.ROI_Abbrev(lutRow);
ROI_EN = LUT.ROI_EN(lutRow);
ROI_CN = LUT.ROI_CN(lutRow);
ROI_TotalVox = LUT.ROI_TotalVox(lutRow);
ROI_SigVox = 0;
ROI_SigRatio = 0;
ROI_SigPercent = 0;
ActivationDegree = 0;
MeanStat_SigVox = NaN;
MedianStat_SigVox = NaN;
MaxStat_SigVox = NaN;
WeightedMeanStat = NaN;
Peak_MNI_X = NaN;
Peak_MNI_Y = NaN;
Peak_MNI_Z = NaN;
Peak_P_unc = NaN;
Peak_Q_FDR = NaN;
T = table(ContrastSet, ContrastName, AnalysisType, StatType, ...
    Merged_ROI_ID, ROI_Abbrev, ROI_EN, ROI_CN, ROI_TotalVox, ...
    ROI_SigVox, ROI_SigRatio, ROI_SigPercent, ActivationDegree, ...
    MeanStat_SigVox, MedianStat_SigVox, MaxStat_SigVox, WeightedMeanStat, ...
    Peak_MNI_X, Peak_MNI_Y, Peak_MNI_Z, Peak_P_unc, Peak_Q_FDR);
end

function Tvox = build_sigvoxel_table(sigIdx, Ystat, pMap, qMap, Yroi, Vstat, LUT, xcon, analysisType, CFG)
n = numel(sigIdx);

if n == 0
    Tvox = empty_voxel_table();
    return;
end

[I,J,K] = ind2sub(size(Ystat), sigIdx);
XYZ = Vstat.mat * [I(:)'; J(:)'; K(:)'; ones(1,n)];
roiID = double(Yroi(sigIdx));
statVal = Ystat(sigIdx);
pVal = pMap(sigIdx);
qVal = qMap(sigIdx);

ContrastSet = repmat({parse_contrastset(xcon.name)}, n, 1);
ContrastName = repmat({xcon.name}, n, 1);
AnalysisType = repmat({analysisType}, n, 1);
StatType = repmat({xcon.STAT}, n, 1);

ROI_ID = roiID(:);
ROI_Abbrev = cell(n,1);
ROI_EN = cell(n,1);
ROI_CN = cell(n,1);
ROI_TotalVox = zeros(n,1);

for ii = 1:n
    lr = find(LUT.ROI_ID == ROI_ID(ii), 1, 'first');
    if isempty(lr)
        ROI_Abbrev{ii} = sprintf('ROI%d', ROI_ID(ii));
        ROI_EN{ii} = ROI_Abbrev{ii};
        ROI_CN{ii} = ROI_Abbrev{ii};
        ROI_TotalVox(ii) = sum(Yroi(:) == ROI_ID(ii));
    else
        ROI_Abbrev{ii} = LUT.ROI_Abbrev{lr};
        ROI_EN{ii} = LUT.ROI_EN{lr};
        ROI_CN{ii} = LUT.ROI_CN{lr};
        ROI_TotalVox(ii) = LUT.ROI_TotalVox(lr);
    end
end

Tvox = table(ContrastSet, ContrastName, AnalysisType, StatType, ...
    statVal(:), pVal(:), qVal(:), ...
    I(:), J(:), K(:), XYZ(1,:)', XYZ(2,:)', XYZ(3,:)', ...
    ROI_ID, ROI_Abbrev, ROI_EN, ROI_CN, ROI_TotalVox, ...
    'VariableNames', { ...
    'ContrastSet','ContrastName','AnalysisType','StatType', ...
    'Stat','P_unc','Q_FDR', ...
    'Voxel_I','Voxel_J','Voxel_K','MNI_X','MNI_Y','MNI_Z', ...
    'Merged_ROI_ID','ROI_Abbrev','ROI_EN','ROI_CN','ROI_TotalVox' ...
    });
end

function T = empty_voxel_table()
T = table(cell(0,1), cell(0,1), cell(0,1), cell(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), cell(0,1), cell(0,1), cell(0,1), zeros(0,1), ...
    'VariableNames', { ...
    'ContrastSet','ContrastName','AnalysisType','StatType', ...
    'Stat','P_unc','Q_FDR', ...
    'Voxel_I','Voxel_J','Voxel_K','MNI_X','MNI_Y','MNI_Z', ...
    'Merged_ROI_ID','ROI_Abbrev','ROI_EN','ROI_CN','ROI_TotalVox' ...
    });
end

function Troi = build_roi_summary_table(sigMask, Ystat, pMap, qMap, Yroi, Vstat, LUT, xcon, analysisType, CFG)
nROI = numel(LUT.ROI_ID);
ContrastSet = repmat({parse_contrastset(xcon.name)}, nROI, 1);
ContrastName = repmat({xcon.name}, nROI, 1);
AnalysisType = repmat({analysisType}, nROI, 1);
StatType = repmat({xcon.STAT}, nROI, 1);

Merged_ROI_ID = LUT.ROI_ID(:);
ROI_Abbrev = LUT.ROI_Abbrev(:);
ROI_EN = LUT.ROI_EN(:);
ROI_CN = LUT.ROI_CN(:);
ROI_TotalVox = LUT.ROI_TotalVox(:);

ROI_SigVox = zeros(nROI,1);
ROI_SigRatio = zeros(nROI,1);
ROI_SigPercent = zeros(nROI,1);
MeanStat_SigVox = nan(nROI,1);
MedianStat_SigVox = nan(nROI,1);
MaxStat_SigVox = nan(nROI,1);
Peak_MNI_X = nan(nROI,1);
Peak_MNI_Y = nan(nROI,1);
Peak_MNI_Z = nan(nROI,1);
Peak_Q_FDR = nan(nROI,1);
Peak_P_unc = nan(nROI,1);
ActivationDegree = zeros(nROI,1);
WeightedMeanStat = nan(nROI,1);

for r = 1:nROI
    id = LUT.ROI_ID(r);
    thisRoi = (Yroi == id) & isfinite(Ystat);
    sig = thisRoi & sigMask;

    ROI_SigVox(r) = sum(sig(:));
    if ROI_TotalVox(r) > 0
        ROI_SigRatio(r) = ROI_SigVox(r) / ROI_TotalVox(r);
    else
        ROI_SigRatio(r) = 0;
    end
    ROI_SigPercent(r) = 100 * ROI_SigRatio(r);
    ActivationDegree(r) = ROI_SigRatio(r);

    if ROI_SigVox(r) > 0
        vals = Ystat(sig);
        qs = qMap(sig);
        ps = pMap(sig);

        MeanStat_SigVox(r) = mean(vals, 'omitnan');
        MedianStat_SigVox(r) = median(vals, 'omitnan');

        [mx, imx] = max(vals);
        MaxStat_SigVox(r) = mx;

        sigLinear = find(sig);
        peakIdx = sigLinear(imx);
        [ii,jj,kk] = ind2sub(size(Ystat), peakIdx);
        xyz = Vstat.mat * [ii; jj; kk; 1];

        Peak_MNI_X(r) = xyz(1);
        Peak_MNI_Y(r) = xyz(2);
        Peak_MNI_Z(r) = xyz(3);
        Peak_Q_FDR(r) = qs(imx);
        Peak_P_unc(r) = ps(imx);
        WeightedMeanStat(r) = ROI_SigRatio(r) * MeanStat_SigVox(r);
    end
end

Troi = table(ContrastSet, ContrastName, AnalysisType, StatType, ...
    Merged_ROI_ID, ROI_Abbrev, ROI_EN, ROI_CN, ROI_TotalVox, ...
    ROI_SigVox, ROI_SigRatio, ROI_SigPercent, ActivationDegree, ...
    MeanStat_SigVox, MedianStat_SigVox, MaxStat_SigVox, WeightedMeanStat, ...
    Peak_MNI_X, Peak_MNI_Y, Peak_MNI_Z, Peak_P_unc, Peak_Q_FDR);
end

function preflight_model_inputs(scanRows, CFG, contrastName, outDir)
fprintf('[PREFLIGHT] Checking model input images and ROI-mask overlap: %s\n', contrastName);

assert(exist(CFG.UnionMaskNii,'file') == 2, 'Missing union mask: %s', CFG.UnionMaskNii);
Vm = spm_vol(CFG.UnionMaskNii);
M = spm_read_vols(Vm) > 0;
roiVox = nnz(M);
if roiVox == 0
    error('Union ROI mask contains 0 voxels: %s', CFG.UnionMaskNii);
end

n = numel(scanRows);
SubjectImage = strings(n,1);
SameDim = false(n,1);
SameMat = false(n,1);
FiniteVoxInROI = zeros(n,1);
NonzeroVoxInROI = zeros(n,1);
MeanAbsInROI = zeros(n,1);
MinInROI = zeros(n,1);
MaxInROI = zeros(n,1);

for ii = 1:n
    f = char(scanRows{ii});
    f = regexprep(f, ',\d+$', '');
    SubjectImage(ii) = string(f);

    if exist(f,'file') ~= 2
        error('Input image missing before SPM estimation: %s', f);
    end

    V = spm_vol(f);
    SameDim(ii) = isequal(V.dim, Vm.dim);
    SameMat(ii) = max(abs(V.mat(:) - Vm.mat(:))) < 1e-4;

    if ~SameDim(ii) || ~SameMat(ii)
        error(['Input image grid does not match ROI mask.\n' ...
               'Image: %s\nMask : %s\n' ...
               'Set RunStep2_CleanReslice=true and SkipExistingAtlasConImages=false to regenerate atlas-space images.'], ...
               f, CFG.UnionMaskNii);
    end

    Y = spm_read_vols(V);
    vals = Y(M);
    finiteVals = vals(isfinite(vals));
    FiniteVoxInROI(ii) = numel(finiteVals);
    NonzeroVoxInROI(ii) = nnz(isfinite(vals) & vals ~= 0);

    if isempty(finiteVals)
        MeanAbsInROI(ii) = NaN;
        MinInROI(ii) = NaN;
        MaxInROI(ii) = NaN;
    else
        MeanAbsInROI(ii) = mean(abs(finiteVals));
        MinInROI(ii) = min(finiteVals);
        MaxInROI(ii) = max(finiteVals);
    end
end

QC = table(SubjectImage, SameDim, SameMat, FiniteVoxInROI, NonzeroVoxInROI, ...
    MeanAbsInROI, MinInROI, MaxInROI);
qcPath = fullfile(outDir, sprintf('Preflight_InputQC_%s.csv', safe_filename(contrastName)));
write_table_utf8bom(QC, qcPath);
fprintf('[WRITE] %s\n', qcPath);
fprintf('[PREFLIGHT] ROI voxels=%d | min nonzero voxels inside ROI across images=%d | median=%g\n', ...
    roiVox, min(NonzeroVoxInROI), median(NonzeroVoxInROI));

if min(FiniteVoxInROI) == 0
    error(['At least one input image has 0 finite voxels inside ROI mask.\n' ...
           'Check Preflight_InputQC CSV. This usually means atlas/resliced images are not in the same valid space.']);
end

if max(NonzeroVoxInROI) < CFG.MinNonzeroVoxInROIMask
    error(['All input images have too few nonzero voxels inside ROI mask.\n' ...
           'This strongly suggests ROI mask and con images do not overlap, or the atlas-space contrast images contain no usable values.\n' ...
           'Contrast images and the ROI atlas must have valid spatial overlap.']);
end

if min(NonzeroVoxInROI) < CFG.MinNonzeroVoxInROIMask
    warning(['Some input images have fewer than %d nonzero voxels inside ROI mask.\n' ...
             'SPM may still estimate the model, but inspect Preflight_InputQC before interpreting results.'], ...
             CFG.MinNonzeroVoxInROIMask);
end
end

function CFG = prepare_functional_roi_atlas(CFG)
% Use the first participant contrast only as a geometry reference. Integer
% labels are sampled with nearest-neighbour interpolation in world space.
contrasts = dir(fullfile(CFG.FirstLevelRoot,'**','con_*.nii'));
assert(~isempty(contrasts),'No first-level contrast images found under %s.',CFG.FirstLevelRoot);
paths = string(fullfile({contrasts.folder},{contrasts.name}));
paths = sort(paths);
referencePath = char(paths(1));
Vref = spm_vol(referencePath);
Vsource = spm_vol(CFG.AtlasNii);
assert(numel(Vref)==1 && numel(Vsource)==1,'The atlas and reference contrast must be 3D volumes.');

labels = zeros(Vref.dim,'uint16');
[x,y] = ndgrid(1:Vref.dim(1),1:Vref.dim(2));
for z = 1:Vref.dim(3)
    referenceVoxels = [x(:)';y(:)';repmat(z,1,numel(x));ones(1,numel(x))];
    sourceVoxels = Vsource.mat \ (Vref.mat*referenceVoxels);
    sampled = spm_sample_vol(Vsource,sourceVoxels(1,:),sourceVoxels(2,:),sourceVoxels(3,:),0);
    sampled(~isfinite(sampled)) = 0;
    labels(:,:,z) = reshape(uint16(round(sampled)),size(x));
end

CFG.AtlasDir = fullfile(CFG.OutRoot,'Atlas');
mkdir_if_needed(CFG.AtlasDir);
CFG.AtlasNii = fullfile(CFG.AtlasDir,'HCPex_Merged17_functional.nii');
CFG.UnionMaskNii = fullfile(CFG.AtlasDir,'HCPex_Merged17_unionmask_functional.nii');
CFG.RoiCountCsv = fullfile(CFG.AtlasDir,'HCPex_Merged17_roi_voxelcount_functional.csv');
Vout = Vref;
Vout.fname = CFG.AtlasNii;
Vout.dt = [spm_type('uint16') 0];
Vout.pinfo = [1;0;0];
Vout.descrip = 'Bilateral HCPex ROIs in functional space: nearest neighbour';
spm_write_vol(Vout,double(labels));
Vmask = Vout;
Vmask.fname = CFG.UnionMaskNii;
Vmask.dt = [spm_type('uint8') 0];
Vmask.descrip = 'HCPex ROI union in functional space';
spm_write_vol(Vmask,double(labels>0));

LUT = load_roi_lookup(CFG);
LUT.ReferenceImage = repmat({referencePath},height(LUT),1);
write_table_utf8bom(LUT,CFG.RoiCountCsv);
fprintf('[ATLAS] Functional-space labels: %s\n',CFG.AtlasNii);
end

function ensure_union_mask(CFG)
if exist(CFG.UnionMaskNii, 'file') == 2
    return;
end
V = spm_vol(CFG.AtlasNii);
Y = spm_read_vols(V);
M = uint8(round(Y) > 0);
Vo = V;
Vo.fname = CFG.UnionMaskNii;
Vo.dt = [2 0];
spm_write_vol(Vo, M);
fprintf('[WRITE] union mask: %s\n', CFG.UnionMaskNii);
end

function LUT = load_roi_lookup(CFG)
T = readtable_safe(CFG.LookupCsv);
colID = pick_col(T, {'Merged_ROI_ID','MergedROI_ID','ROI_ID','MergedID','Merged_ROI'});
colAbbrev = pick_col_optional(T, {'Merged_ROI_Abbrev','MergedROI_Abbrev','ROI_Abbrev','Abbrev'});
colEN = pick_col_optional(T, {'Merged_ROI_EN','MergedROI_EN','ROI_EN','ROI_Name_EN','EN'});
colCN = pick_col_optional(T, {'Merged_ROI_CN','MergedROI_CN','ROI_CN','ROI_Name_CN','CN'});

ids = double(T.(colID));
uID = unique(ids(isfinite(ids) & ids > 0), 'stable');

ROI_ID = zeros(numel(uID),1);
ROI_Abbrev = cell(numel(uID),1);
ROI_EN = cell(numel(uID),1);
ROI_CN = cell(numel(uID),1);
ROI_TotalVox = zeros(numel(uID),1);

V = spm_vol(CFG.AtlasNii);
Y = round(spm_read_vols(V));

for i = 1:numel(uID)
    id = uID(i);
    row = find(ids == id, 1, 'first');
    ROI_ID(i) = id;

    if ~isempty(colAbbrev), ROI_Abbrev{i} = to_char(T.(colAbbrev)(row)); else, ROI_Abbrev{i} = sprintf('ROI%02d', id); end
    if ~isempty(colEN), ROI_EN{i} = to_char(T.(colEN)(row)); else, ROI_EN{i} = ROI_Abbrev{i}; end
    if ~isempty(colCN), ROI_CN{i} = to_char(T.(colCN)(row)); else, ROI_CN{i} = ROI_Abbrev{i}; end

    ROI_TotalVox(i) = sum(Y(:) == id);
end

LUT = table(ROI_ID, ROI_Abbrev, ROI_EN, ROI_CN, ROI_TotalVox);
end

function taskName = task_from_group_run(groupName, runName)
g = sscanf(groupName, 'Group%d');
r = sscanf(runName, 'Run%d');

tasks = {'Safety','Vitality','Beauty','Imageability'};

map = [1 2 3 4; ...
       2 1 4 3; ...
       3 4 1 2; ...
       4 3 2 1];

taskName = tasks{map(g,r)};
end

function idx = task_index(taskName, CFG)
idx = find(strcmp(CFG.Tasks, taskName), 1, 'first');
if isempty(idx)
    error('Unknown task: %s', taskName);
end
end

function mkdir_if_needed(p)
if exist(p,'dir') ~= 7
    mkdir(p);
end
end

function T = readtable_safe(csvPath)
try
    opts = detectImportOptions(csvPath, 'Encoding','UTF-8');
    opts.VariableNamingRule = 'preserve';
    T = readtable(csvPath, opts);
catch
    T = readtable(csvPath, 'FileType','text', 'Encoding','UTF-8', 'VariableNamingRule','preserve');
end
T.Properties.VariableNames = matlab.lang.makeValidName(T.Properties.VariableNames, 'ReplacementStyle','delete');
end

function col = pick_col(T, candidates)
vars = T.Properties.VariableNames;
for i = 1:numel(candidates)
    c = matlab.lang.makeValidName(candidates{i}, 'ReplacementStyle','delete');
    hit = strcmp(vars, c);
    if any(hit)
        col = vars{find(hit,1,'first')};
        return;
    end
end
error('Cannot find required column. Tried: %s', strjoin(candidates, ', '));
end

function col = pick_col_optional(T, candidates)
vars = T.Properties.VariableNames;
col = '';
for i = 1:numel(candidates)
    c = matlab.lang.makeValidName(candidates{i}, 'ReplacementStyle','delete');
    hit = strcmp(vars, c);
    if any(hit)
        col = vars{find(hit,1,'first')};
        return;
    end
end
end

function s = to_char(x)
if iscell(x), x = x{1}; end
if isstring(x), s = char(x); return; end
if ischar(x), s = x; return; end
if isnumeric(x), s = num2str(x); return; end
try
    s = char(x);
catch
    s = '';
end
end

function write_table_utf8bom(T, outPath)
outDir = fileparts(outPath);
if ~isempty(outDir), mkdir_if_needed(outDir); end

tmp = [tempname, '.csv'];

try
    writetable(T, tmp, 'FileType','text', 'Encoding','UTF-8');
catch
    writetable(T, tmp, 'Encoding','UTF-8');
end

fidIn = fopen(tmp, 'rb');
if fidIn < 0
    error('Cannot open temporary CSV for reading: %s', tmp);
end
bytes = fread(fidIn, '*uint8');
fclose(fidIn);

bom = uint8([239 187 191]);
if numel(bytes) >= 3 && all(bytes(1:3)' == bom)
    bytes = bytes(4:end);
end

fidOut = fopen(outPath, 'wb');
if fidOut < 0
    error('Cannot open CSV for writing: %s', outPath);
end
fwrite(fidOut, bom, 'uint8');  % UTF-8 BOM for Excel Chinese compatibility
fwrite(fidOut, bytes, 'uint8');
fclose(fidOut);

try, delete(tmp); catch, end
end

function q = bh_fdr_q(p)
p = p(:);
m = numel(p);
[ps, ord] = sort(p, 'ascend');
qSorted = ps .* m ./ (1:m)';
qSorted = flipud(cummin(flipud(qSorted)));
qSorted(qSorted > 1) = 1;
q = nan(size(p));
q(ord) = qSorted;
end

function s = safe_name(s)
s = safe_filename(s);
end

function s = safe_filename(s)
s = char(s);
s = regexprep(s, '[<>:"/\\|?*\s]+', '_');
s = regexprep(s, '_+', '_');
end

function cs = parse_contrastset(cname)
cname = char(cname);
parts = split(string(cname), "_");
if numel(parts) >= 3 && startsWith(cname, 'Stim_gt_Fixation')
    cs = 'Stim_gt_Fixation';
elseif startsWith(cname, 'RespBin_gt_0')
    cs = 'RespBin_gt_0';
elseif startsWith(cname, 'RespBin_lt_0')
    cs = 'RespBin_lt_0';
else
    cs = char(parts(1));
end
end

function atype = classify_contrast_name(cname)
cname = char(cname);
if contains(cname, '_OneSample_')
    atype = 'OneSample';
elseif contains(cname, '_TaskMainEffect_F')
    atype = 'TaskMainEffect';
elseif contains(cname, '_TaskXGroup_Interaction_F')
    atype = 'TaskXGroupInteraction';
elseif contains(cname, '_neq_F') || contains(cname, '_vs_') && endsWith(cname, '_F')
    atype = 'PairwiseEqualityF';
elseif contains(cname, '_gt_')
    atype = 'Pairwise';
else
    atype = 'Unknown';
end
end
