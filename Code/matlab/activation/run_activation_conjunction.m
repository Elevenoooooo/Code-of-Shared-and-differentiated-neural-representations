function run_activation_conjunction(cfg)
% RUN_ACTIVATION_CONJUNCTION Intersect the four positive task activation maps.
% Each input is voxelwise-FDR significant at q<0.05. The display statistic is
% the minimum T across tasks within the strict intersection of their masks.
% The final ROI set excludes the amygdala and entorhinal cortex.
% Overlap bit codes: Safety=1, Vitality=2, Beauty=4, Imageability=8.

if nargin < 1, cfg = study_config(); end
required = {'spm_dir','activation_root','merged_atlas_dir','hcpex_atlas','hcpex_labels'};
for k = 1:numel(required)
    assert(isfield(cfg,required{k}) && ~isempty(cfg.(required{k})), 'Missing cfg.%s', required{k});
end
addpath(cfg.spm_dir);

spm('Defaults','fMRI');
spm_jobman('initcfg');

CFG.OutRoot = cfg.activation_root;

CFG.AtlasDir = cfg.merged_atlas_dir;
CFG.AtlasNii = cfg.hcpex_atlas;
CFG.LookupCsv = cfg.hcpex_labels;
functionalAtlas = fullfile(cfg.activation_root,'Atlas','HCPex_Merged17_functional.nii');
if isfile(functionalAtlas)
    CFG.AtlasNii = functionalAtlas;
    functionalCounts = fullfile(cfg.activation_root,'Atlas','HCPex_Merged17_roi_voxelcount_functional.csv');
    if isfile(functionalCounts), CFG.LookupCsv = functionalCounts; end
end

CFG.ContrastSets = {'Stim_gt_Fixation'};
CFG.Tasks = {'Safety','Vitality','Beauty','Imageability'};

CFG.FDR_q = 0.05;
CFG.MapDirName = 'FDR_Maps';
CFG.OutDirName = 'FourTask_ActivationConjunction';

CFG.Use15ROI_ExcludeAmygEntorhinal = true;

CFG.WriteGe3TasksMask = false;

if isfield(cfg,'conjunction_options')
    user = cfg.conjunction_options; names = fieldnames(user);
    for k = 1:numel(names)
        assert(isfield(CFG,names{k}), 'Unknown conjunction_options option: %s',names{k});
        CFG.(names{k}) = user.(names{k});
    end
end
assert(exist(CFG.OutRoot,'dir') == 7, 'Missing OutRoot: %s', CFG.OutRoot);
assert(exist(CFG.AtlasNii,'file') == 2, 'Missing AtlasNii: %s', CFG.AtlasNii);

mapDir = fullfile(CFG.OutRoot, CFG.MapDirName);
assert(exist(mapDir,'dir') == 7, 'Missing FDR_Maps dir: %s', mapDir);

outRoot = fullfile(CFG.OutRoot, CFG.OutDirName);
mkdir_if_needed(outRoot);

Vroi = spm_vol(CFG.AtlasNii);
Yroi = round(spm_read_vols(Vroi));
roiMask0 = Yroi > 0 & isfinite(Yroi);

LUT = load_roi_lookup_for_conjunction(CFG, Yroi);
keepROI = get_keep_roi_flags(LUT, CFG);

keepIDs = LUT.ROI_ID(keepROI);
roiMask = roiMask0 & ismember(Yroi, keepIDs);

fprintf('[ROI] total ROI labels in lookup = %d\n', numel(LUT.ROI_ID));
fprintf('[ROI] used ROI labels = %d\n', numel(keepIDs));
fprintf('[ROI] used ROI voxels = %d\n', nnz(roiMask));

for ci = 1:numel(CFG.ContrastSets)
    csName = CFG.ContrastSets{ci};
    fprintf('\n============================================================\n');
    fprintf('[CONJ] Four-task activation conjunction: %s\n', csName);
    fprintf('============================================================\n');

    outDir = fullfile(outRoot, csName);
    mkdir_if_needed(outDir);

    nT = numel(CFG.Tasks);
    sigMasks = cell(nT,1);
    statMaps = cell(nT,1);
    Vref = [];

    for ti = 1:nT
        taskName = CFG.Tasks{ti};
        [sigPath, statPath] = find_task_activation_maps(mapDir, csName, taskName, CFG);

        fprintf('[READ] %s\n', sigPath);
        Vsig = spm_vol(sigPath);
        Ysig = spm_read_vols(Vsig) > 0;

        fprintf('[READ] %s\n', statPath);
        Vstat = spm_vol(statPath);
        Ystat = spm_read_vols(Vstat);

        if isempty(Vref)
            Vref = Vstat;
            assert_same_grid(Vroi, Vref, 'Atlas and first stat map');
        else
            assert_same_grid(Vref, Vstat, sprintf('Reference and stat map: %s', taskName));
        end
        assert_same_grid(Vref, Vsig, sprintf('Reference and sig mask: %s', taskName));

        sigMasks{ti} = Ysig & roiMask;
        statMaps{ti} = Ystat;
    end

    overlapCount = zeros(size(Yroi), 'uint8');
    labelMap = zeros(size(Yroi), 'uint8');

    for ti = 1:nT
        m = sigMasks{ti} & roiMask;
        overlapCount(m) = overlapCount(m) + uint8(1);
        labelMap(m) = labelMap(m) + uint8(2^(ti-1));
    end

    commonMask = roiMask & (overlapCount == nT);

    minT = statMaps{1};
    for ti = 2:nT
        minT = min(minT, statMaps{ti});
    end

    minTCommon = zeros(size(minT), 'single');
    minTCommon(commonMask) = single(minT(commonMask));

    ge3Mask = roiMask & (overlapCount >= 3);

    fprintf('[RESULT] common 4-task voxels = %d\n', nnz(commonMask));
    fprintf('[RESULT] >=3-task overlap voxels = %d\n', nnz(ge3Mask));

    prefix = safe_filename(sprintf('%s_FourTask', csName));

    commonMaskPath = fullfile(outDir, sprintf('%s_Common4Task_Activation_mask.nii', prefix));
    minTPath       = fullfile(outDir, sprintf('%s_Common4Task_Activation_minT.nii', prefix));
    countPath      = fullfile(outDir, sprintf('%s_TaskOverlap_CountMap.nii', prefix));
    labelPath      = fullfile(outDir, sprintf('%s_TaskOverlap_LabelMap_BitCode.nii', prefix));
    ge3Path        = fullfile(outDir, sprintf('%s_TaskOverlap_Ge3Tasks_mask.nii', prefix));

    write_nifti_like(Vref, uint8(commonMask), commonMaskPath, [2 0]);
    write_nifti_like(Vref, minTCommon, minTPath, [16 0]);
    write_nifti_like(Vref, uint8(overlapCount), countPath, [2 0]);
    write_nifti_like(Vref, uint8(labelMap), labelPath, [2 0]);

    if CFG.WriteGe3TasksMask
        write_nifti_like(Vref, uint8(ge3Mask), ge3Path, [2 0]);
    end

    fprintf('[WRITE] %s\n', commonMaskPath);
    fprintf('[WRITE] %s\n', minTPath);
    fprintf('[WRITE] %s\n', countPath);
    fprintf('[WRITE] %s\n', labelPath);
    if CFG.WriteGe3TasksMask
        fprintf('[WRITE] %s\n', ge3Path);
    end

    Troi = build_common4_roi_summary(commonMask, overlapCount, minT, Yroi, Vref, LUT, keepROI, csName, CFG);
    outRoiCsv = fullfile(outDir, sprintf('%s_Common4Task_ROISummary.csv', prefix));
    write_table_utf8bom(Troi, outRoiCsv);
    fprintf('[WRITE] %s\n', outRoiCsv);

    TcombGlobal = build_combination_global_summary(labelMap, roiMask, CFG, csName);
    outCombGlobal = fullfile(outDir, sprintf('%s_TaskCombination_GlobalSummary.csv', prefix));
    write_table_utf8bom(TcombGlobal, outCombGlobal);
    fprintf('[WRITE] %s\n', outCombGlobal);

    TcombROI = build_combination_byroi_summary(labelMap, Yroi, LUT, keepROI, CFG, csName);
    outCombRoi = fullfile(outDir, sprintf('%s_TaskCombination_ByROISummary.csv', prefix));
    write_table_utf8bom(TcombROI, outCombRoi);
    fprintf('[WRITE] %s\n', outCombRoi);
end

fprintf('\n[ALL DONE] Four-task activation conjunction finished.\n');

end

function [sigPath, statPath] = find_task_activation_maps(mapDir, csName, taskName, CFG)

qtag = sprintf('q%03d', round(CFG.FDR_q * 1000));
base = safe_filename(sprintf('%s_OneSample_%s_gt_0', csName, taskName));

sigPath = fullfile(mapDir, sprintf('%s_sigmask_%s.nii', base, qtag));
if exist(sigPath, 'file') ~= 2
    error(['Missing one-sample FDR sigmask:\n  %s\n\n' ...
           'Generate task one-sample FDR outputs with run_activation_group_model:\n' ...
           '  CFG.FDR_IncludeTaskOneSample = true;\n' ...
           '  CFG.RunStep3A_TaskOneSample = true;\n' ...
           '  CFG.RunStep4_FDROutputs = true;\n'], sigPath);
end

statPath = fullfile(mapDir, sprintf('%s_statmap.nii', base));
if exist(statPath, 'file') ~= 2
    statPath2 = fullfile(mapDir, sprintf('%s_statmap_%s.nii', base, qtag));
    if exist(statPath2, 'file') == 2
        statPath = statPath2;
    else
        error('Missing stat map for %s | %s:\n  %s\n  %s', csName, taskName, statPath, statPath2);
    end
end

end

function Troi = build_common4_roi_summary(commonMask, overlapCount, minT, Yroi, Vref, LUT, keepROI, csName, CFG)

ids = LUT.ROI_ID(keepROI);
nROI = numel(ids);

ContrastSet = repmat({csName}, nROI, 1);
ConjunctionName = repmat({'Safety_AND_Vitality_AND_Beauty_AND_Imageability_gt0_qFDR05'}, nROI, 1);
Tasks = repmat({strjoin(CFG.Tasks, '|')}, nROI, 1);
NTasks = repmat(numel(CFG.Tasks), nROI, 1);

Merged_ROI_ID = ids(:);
ROI_Abbrev = LUT.ROI_Abbrev(keepROI);
ROI_EN = LUT.ROI_EN(keepROI);
ROI_CN = LUT.ROI_CN(keepROI);

ROI_TotalVox_Used = zeros(nROI,1);
Common4_Vox = zeros(nROI,1);
Common4_Percent = zeros(nROI,1);

Overlap1_Vox = zeros(nROI,1);
Overlap2_Vox = zeros(nROI,1);
Overlap3_Vox = zeros(nROI,1);
Overlap4_Vox = zeros(nROI,1);

Peak_MinT = nan(nROI,1);
Peak_MNI_X = nan(nROI,1);
Peak_MNI_Y = nan(nROI,1);
Peak_MNI_Z = nan(nROI,1);

Mean_MinT_Common = nan(nROI,1);
Median_MinT_Common = nan(nROI,1);

for r = 1:nROI
    id = ids(r);
    roi = (Yroi == id);
    ROI_TotalVox_Used(r) = nnz(roi);

    Common4_Vox(r) = nnz(roi & commonMask);
    if ROI_TotalVox_Used(r) > 0
        Common4_Percent(r) = 100 * Common4_Vox(r) / ROI_TotalVox_Used(r);
    end

    Overlap1_Vox(r) = nnz(roi & overlapCount == 1);
    Overlap2_Vox(r) = nnz(roi & overlapCount == 2);
    Overlap3_Vox(r) = nnz(roi & overlapCount == 3);
    Overlap4_Vox(r) = nnz(roi & overlapCount == 4);

    thisCommon = roi & commonMask;
    if any(thisCommon(:))
        vals = minT(thisCommon);
        Mean_MinT_Common(r) = mean(vals, 'omitnan');
        Median_MinT_Common(r) = median(vals, 'omitnan');

        [mx, imx] = max(vals);
        Peak_MinT(r) = mx;

        idx = find(thisCommon);
        peakIdx = idx(imx);
        [ii,jj,kk] = ind2sub(size(Yroi), peakIdx);
        xyz = Vref.mat * [ii; jj; kk; 1];

        Peak_MNI_X(r) = xyz(1);
        Peak_MNI_Y(r) = xyz(2);
        Peak_MNI_Z(r) = xyz(3);
    end
end

Troi = table( ...
    ContrastSet, ConjunctionName, Tasks, NTasks, ...
    Merged_ROI_ID, ROI_Abbrev, ROI_EN, ROI_CN, ROI_TotalVox_Used, ...
    Common4_Vox, Common4_Percent, ...
    Overlap1_Vox, Overlap2_Vox, Overlap3_Vox, Overlap4_Vox, ...
    Mean_MinT_Common, Median_MinT_Common, Peak_MinT, ...
    Peak_MNI_X, Peak_MNI_Y, Peak_MNI_Z);

end

function T = build_combination_global_summary(labelMap, roiMask, CFG, csName)

nT = numel(CFG.Tasks);
codes = (1:(2^nT - 1))';

ContrastSet = repmat({csName}, numel(codes), 1);
Code = codes;
Combination = cell(numel(codes),1);
NTasks = zeros(numel(codes),1);
VoxelCount = zeros(numel(codes),1);
PercentOfROIMask = zeros(numel(codes),1);

totalVox = nnz(roiMask);

for i = 1:numel(codes)
    code = codes(i);
    taskFlags = false(1,nT);
    for ti = 1:nT
        taskFlags(ti) = bitget(code, ti) == 1;
    end

    Combination{i} = strjoin(CFG.Tasks(taskFlags), '_AND_');
    NTasks(i) = sum(taskFlags);

    m = roiMask & labelMap == code;
    VoxelCount(i) = nnz(m);
    if totalVox > 0
        PercentOfROIMask(i) = 100 * VoxelCount(i) / totalVox;
    end
end

T = table(ContrastSet, Code, Combination, NTasks, VoxelCount, PercentOfROIMask);

end

function T = build_combination_byroi_summary(labelMap, Yroi, LUT, keepROI, CFG, csName)

ids = LUT.ROI_ID(keepROI);
nROI = numel(ids);
nT = numel(CFG.Tasks);
codes = 1:(2^nT - 1);

rows = {};

for r = 1:nROI
    id = ids(r);
    roi = Yroi == id;
    roiTotal = nnz(roi);

    for ci = 1:numel(codes)
        code = codes(ci);

        taskFlags = false(1,nT);
        for ti = 1:nT
            taskFlags(ti) = bitget(code, ti) == 1;
        end

        comboName = strjoin(CFG.Tasks(taskFlags), '_AND_');
        nTasks = sum(taskFlags);

        vox = nnz(roi & labelMap == code);
        pct = 0;
        if roiTotal > 0
            pct = 100 * vox / roiTotal;
        end

        rows(end+1,:) = { ...
            csName, id, LUT.ROI_Abbrev{keepROI_index(keepROI, r)}, ...
            LUT.ROI_EN{keepROI_index(keepROI, r)}, LUT.ROI_CN{keepROI_index(keepROI, r)}, ...
            roiTotal, code, comboName, nTasks, vox, pct}; %#ok<AGROW>
    end
end

T = cell2table(rows, 'VariableNames', { ...
    'ContrastSet','Merged_ROI_ID','ROI_Abbrev','ROI_EN','ROI_CN','ROI_TotalVox_Used', ...
    'Code','Combination','NTasks','VoxelCount','PercentOfROI'});

end

function idx = keepROI_index(keepROI, r)
allIdx = find(keepROI);
idx = allIdx(r);
end

function LUT = load_roi_lookup_for_conjunction(CFG, Yroi)

if exist(CFG.LookupCsv, 'file') ~= 2
    warning('LookupCsv not found. ROI names will be generated from atlas labels.');
    ids = unique(Yroi(Yroi > 0));
    LUT.ROI_ID = double(ids(:));
    LUT.ROI_Abbrev = cellstr("ROI" + string(ids(:)));
    LUT.ROI_EN = LUT.ROI_Abbrev;
    LUT.ROI_CN = LUT.ROI_Abbrev;
    LUT.ROI_TotalVox = arrayfun(@(x) nnz(Yroi == x), ids(:));
    return;
end

T = readtable(CFG.LookupCsv, 'TextType','string', 'VariableNamingRule','preserve');
vars = string(T.Properties.VariableNames);

idCol = pick_col(vars, {'Merged_ROI_ID','ROI_ID','ROIID','ID','Label','Value'});
if isempty(idCol)
    error('Cannot identify ROI ID column in lookup CSV: %s', CFG.LookupCsv);
end

abbrCol = pick_col(vars, {'ROI_Abbrev','Merged_ROI_Abbrev','Abbrev','ROIAbbrev','ROI'});
enCol    = pick_col(vars, {'ROI_EN','Merged_ROI_EN','ROI_Name_EN','ROIName','Name','ROI_EN_Name'});
cnCol    = pick_col(vars, {'ROI_CN','Merged_ROI_CN','ROI_Name_CN','ChineseName','CN'});
voxCol   = pick_col(vars, {'ROI_TotalVox','Merged_ROI_TotalVox','TotalVox','VoxelCount','Nvox','ROI_Voxels'});

LUT.ROI_ID = double(T.(char(vars(idCol))));

if isempty(abbrCol)
    LUT.ROI_Abbrev = cellstr("ROI" + string(LUT.ROI_ID));
else
    LUT.ROI_Abbrev = cellstr(string(T.(char(vars(abbrCol)))));
end

if isempty(enCol)
    LUT.ROI_EN = LUT.ROI_Abbrev;
else
    LUT.ROI_EN = cellstr(string(T.(char(vars(enCol)))));
end

if isempty(cnCol)
    LUT.ROI_CN = LUT.ROI_Abbrev;
else
    LUT.ROI_CN = cellstr(string(T.(char(vars(cnCol)))));
end

if isempty(voxCol)
    LUT.ROI_TotalVox = arrayfun(@(x) nnz(Yroi == x), LUT.ROI_ID(:));
else
    LUT.ROI_TotalVox = double(T.(char(vars(voxCol))));
end

valid = isfinite(LUT.ROI_ID) & LUT.ROI_ID > 0;
LUT.ROI_ID = LUT.ROI_ID(valid);
LUT.ROI_Abbrev = LUT.ROI_Abbrev(valid);
LUT.ROI_EN = LUT.ROI_EN(valid);
LUT.ROI_CN = LUT.ROI_CN(valid);
LUT.ROI_TotalVox = LUT.ROI_TotalVox(valid);

end

function idx = pick_col(vars, candidates)
idx = [];
nv = normalize_name(vars);
nc = normalize_name(string(candidates));
for i = 1:numel(nc)
    hit = find(nv == nc(i), 1, 'first');
    if ~isempty(hit)
        idx = hit;
        return;
    end
end
end

function s = normalize_name(s)
s = string(s);
s = replace(s, char(65279), "");
s = lower(regexprep(s, '[^A-Za-z0-9]', ''));
end

function keepROI = get_keep_roi_flags(LUT, CFG)

keepROI = true(numel(LUT.ROI_ID), 1);

if ~isfield(CFG, 'Use15ROI_ExcludeAmygEntorhinal') || ~CFG.Use15ROI_ExcludeAmygEntorhinal
    return;
end

abbr = lower(strtrim(string(LUT.ROI_Abbrev(:))));
en   = lower(strtrim(string(LUT.ROI_EN(:))));
cn   = lower(strtrim(string(LUT.ROI_CN(:))));
txt  = abbr + " " + en + " " + cn;

isAmyg = contains(txt, "amyg") | contains(txt, "amygdala") | contains(txt, "杏仁");
isEnt  = contains(txt, "entorhinal") | contains(txt, "内嗅");

isEC = abbr == "ec" | abbr == "erc" | abbr == "entorhinal cortex";

keepROI = ~(isAmyg | isEnt | isEC);

fprintf('[ROI] Excluding Amygdala / Entorhinal cortex: nExcluded = %d\n', nnz(~keepROI));

end

function write_nifti_like(Vref, Y, outPath, dt)
Vout = Vref;
Vout.fname = outPath;
Vout.dt = dt;
spm_write_vol(Vout, Y);
end

function assert_same_grid(Va, Vb, label)
sameDim = isequal(Va.dim, Vb.dim);
sameMat = max(abs(Va.mat(:) - Vb.mat(:))) < 1e-4;
if ~(sameDim && sameMat)
    error('Grid mismatch: %s', label);
end
end

function mkdir_if_needed(d)
if exist(d, 'dir') ~= 7
    mkdir(d);
end
end

function s = safe_filename(s)
s = char(string(s));
s = regexprep(s, '[^\w\-\.\(\)]', '_');
end

function write_table_utf8bom(T, outPath)
writetable(T, outPath, ...
    'FileType','text', ...
    'Delimiter', ',', ...
    'Encoding', 'UTF-8', ...
    'WriteVariableNames', true);
end
