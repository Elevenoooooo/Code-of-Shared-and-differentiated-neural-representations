function lookupOutput = update_atlas_voxel_counts(cfg)
% UPDATE_ATLAS_VOXEL_COUNTS Count voxels for each merged ROI in atlas space.
% Reads cfg.hcpex_atlas and cfg.hcpex_labels, writes a derivative lookup with
% Merged_ROI_TotalVox, and returns its path for subsequent analysis stages.

if nargin < 1, cfg = study_config(); end
required = {'spm_dir','merged_atlas_dir','hcpex_atlas','hcpex_labels'};
for k = 1:numel(required)
    assert(isfield(cfg,required{k}) && ~isempty(cfg.(required{k})), 'Missing cfg.%s', required{k});
end
addpath(cfg.spm_dir);

spm('Defaults','fMRI');

AtlasDir = cfg.merged_atlas_dir;

AtlasNii = cfg.hcpex_atlas;
LookupCsv = cfg.hcpex_labels;
lookupOutput = fullfile(AtlasDir,'HCPex_MergedROI_UrbanPerception_Refined17_lookup.csv');
if ~isfolder(AtlasDir), mkdir(AtlasDir); end

OutRoiCountCsv = fullfile(AtlasDir, 'HCPex_MergedROI_UrbanPerception_Refined17_roi_voxelcount.csv');

BackupOriginalLookup = true;

assert(exist(AtlasNii,'file') == 2, 'Missing atlas NIfTI: %s', AtlasNii);
assert(exist(LookupCsv,'file') == 2, 'Missing lookup CSV: %s', LookupCsv);

V = spm_vol(AtlasNii);
Y = spm_read_vols(V);
Y = round(Y);
roiIDsInAtlas = unique(Y(:));
roiIDsInAtlas = roiIDsInAtlas(isfinite(roiIDsInAtlas) & roiIDsInAtlas > 0);

fprintf('[INFO] Atlas: %s\n', AtlasNii);
fprintf('[INFO] Found %d non-zero merged ROI labels.\n', numel(roiIDsInAtlas));

T = readtable_safe(LookupCsv);
colID = pick_col(T, {'Merged_ROI_ID','MergedROI_ID','ROI_ID','MergedID','Merged_ROI'});
colAbbrev = pick_col_optional(T, {'Merged_ROI_Abbrev','MergedROI_Abbrev','ROI_Abbrev','Abbrev'});
colEN = pick_col_optional(T, {'Merged_ROI_EN','MergedROI_EN','ROI_EN','ROI_Name_EN','EN'});
colCN = pick_col_optional(T, {'Merged_ROI_CN','MergedROI_CN','ROI_CN','ROI_Name_CN','CN'});

lookupID = double(T.(colID));
mergedIDs = unique(lookupID(isfinite(lookupID) & lookupID > 0), 'stable');

nROI = numel(mergedIDs);
ROI_ID = zeros(nROI,1);
ROI_Abbrev = cell(nROI,1);
ROI_EN = cell(nROI,1);
ROI_CN = cell(nROI,1);
ROI_TotalVox = zeros(nROI,1);

for r = 1:nROI
    id = mergedIDs(r);
    ROI_ID(r) = id;
    ROI_TotalVox(r) = sum(Y(:) == id);

    firstRow = find(lookupID == id, 1, 'first');

    if ~isempty(colAbbrev)
        ROI_Abbrev{r} = to_char(T.(colAbbrev)(firstRow));
    else
        ROI_Abbrev{r} = sprintf('ROI%02d', id);
    end

    if ~isempty(colEN)
        ROI_EN{r} = to_char(T.(colEN)(firstRow));
    else
        ROI_EN{r} = ROI_Abbrev{r};
    end

    if ~isempty(colCN)
        ROI_CN{r} = to_char(T.(colCN)(firstRow));
    else
        ROI_CN{r} = ROI_Abbrev{r};
    end

    fprintf('[ROI] ID=%d | %-12s | vox=%d\n', id, ROI_Abbrev{r}, ROI_TotalVox(r));
end

Tout = table(ROI_ID, ROI_Abbrev, ROI_EN, ROI_CN, ROI_TotalVox);
writetable(Tout, OutRoiCountCsv, 'Encoding','UTF-8');
fprintf('[WRITE] %s\n', OutRoiCountCsv);

if BackupOriginalLookup && isfile(lookupOutput)
    ts = datestr(now, 'yyyymmdd_HHMMSS');
    backupCsv = fullfile(AtlasDir, ['HCPex_MergedROI_UrbanPerception_Refined17_lookup_backup_' ts '.csv']);
    copyfile(lookupOutput, backupCsv);
    fprintf('[BACKUP] %s\n', backupCsv);
end

totalVoxPerLookupRow = zeros(height(T),1);
for i = 1:height(T)
    id = double(T.(colID)(i));
    idx = find(ROI_ID == id, 1, 'first');
    if ~isempty(idx)
        totalVoxPerLookupRow(i) = ROI_TotalVox(idx);
    else
        totalVoxPerLookupRow(i) = 0;
    end
end

T.Merged_ROI_TotalVox = totalVoxPerLookupRow;

writetable(T, lookupOutput, 'Encoding','UTF-8');
fprintf('[WRITE] Derivative lookup with Merged_ROI_TotalVox:\n%s\n', lookupOutput);
fprintf('[DONE]\n');

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
