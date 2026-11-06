function build_merged_hcpex_atlas(cfg)
% BUILD_MERGED_HCPEX_ATLAS Merge HCPex parcels into 17 bilateral regions of interest.
% The workbook supplies explicit Parcel_Map and ROI_Definition tables.
% Parcel assignments are mutually exclusive. Outputs include the integer
% atlas, parcel lookup, ROI voxel counts and a compact metadata MAT file.

if nargin < 1, cfg = study_config(); end
required = {'spm_dir','hcpex_source_atlas','merged_atlas_dir','atlas_crosswalk'};
for k = 1:numel(required)
    assert(isfield(cfg,required{k}) && ~isempty(cfg.(required{k})), 'Missing cfg.%s', required{k});
end
addpath(cfg.spm_dir);

spm('Defaults','fMRI');
spm_jobman('initcfg');

CFG = struct();

CFG.HCPexAtlasNii = cfg.hcpex_source_atlas;
CFG.OutDir = cfg.merged_atlas_dir;
CFG.CrosswalkXlsx = cfg.atlas_crosswalk;

CFG.OutNii       = fullfile(CFG.OutDir, 'HCPex_MergedROI_UrbanPerception_Refined17_Bilateral_2mm.nii');
CFG.OutLookupCsv = fullfile(CFG.OutDir, 'HCPex_MergedROI_UrbanPerception_Refined17_lookup.csv');
CFG.OutVoxelCsv  = fullfile(CFG.OutDir, 'HCPex_MergedROI_UrbanPerception_Refined17_roi_voxelcount.csv');
CFG.OutCacheMAT  = fullfile(CFG.OutDir, 'HCPex_MergedROI_UrbanPerception_Refined17_Cache.mat');
CFG.Overwrite = false; % Set cfg.atlas_options.Overwrite=true to rebuild explicitly.

if isfield(cfg,'atlas_options')
    user = cfg.atlas_options; names = fieldnames(user);
    for k = 1:numel(names)
        assert(isfield(CFG,names{k}), 'Unknown atlas_options option: %s',names{k});
        CFG.(names{k}) = user.(names{k});
    end
end

mkdir_if_needed(CFG.OutDir);

fprintf('\n====================================================\n');
fprintf('[S03] Build HCPex Merged17 bilateral atlas\n');
fprintf('HCPex atlas  : %s\n', CFG.HCPexAtlasNii);
fprintf('Crosswalk    : %s\n', CFG.CrosswalkXlsx);
fprintf('Output atlas : %s\n', CFG.OutNii);
fprintf('====================================================\n');

if ~exist(CFG.HCPexAtlasNii, 'file')
    error('Cannot find HCPex atlas: %s', CFG.HCPexAtlasNii);
end
if ~exist(CFG.CrosswalkXlsx, 'file')
    error('Cannot find crosswalk xlsx: %s', CFG.CrosswalkXlsx);
end
if exist(CFG.OutNii, 'file') && ~CFG.Overwrite
    fprintf('[SKIP] Output atlas already exists and CFG.Overwrite=false.\n');
    return;
end

Tmap = readtable_preserve(CFG.CrosswalkXlsx, 'Sheet', 'Parcel_Map');
Tdef = readtable_preserve(CFG.CrosswalkXlsx, 'Sheet', 'ROI_Definition');

Tmap = normalize_table_columns(Tmap);
Tdef = normalize_table_columns(Tdef);

vn_map = Tmap.Properties.VariableNames;
vn_def = Tdef.Properties.VariableNames;

idCol       = pick_col_robust(vn_map, {'HCPex_ID'});
lrCol       = pick_col_robust(vn_map, {'HCPex_LR'});
parcelCol   = pick_col_robust(vn_map, {'HCPex_Parcel'});
roiIdCol    = pick_col_robust(vn_map, {'Merged_ROI_ID','ROI_ID'});
roiAbbrCol  = pick_col_robust(vn_map, {'Merged_ROI_Abbrev','ROI_Abbrev'});
roiEnCol    = pick_col_robust(vn_map, {'Merged_ROI_EN','ROI_FullName_EN'});
roiCnCol    = pick_col_robust(vn_map, {'Merged_ROI_CN','ROI_CN'});

needCols = {idCol, lrCol, parcelCol, roiIdCol, roiAbbrCol, roiEnCol, roiCnCol};
if any(cellfun(@isempty, needCols))
    error(['Parcel_Map sheet must contain: HCPex_ID, HCPex_LR, HCPex_Parcel, ', ...
           'Merged_ROI_ID, Merged_ROI_Abbrev, Merged_ROI_EN, Merged_ROI_CN.']);
end

defRoiIdCol   = pick_col_robust(vn_def, {'ROI_ID','Merged_ROI_ID'});
defRoiEnCol   = pick_col_robust(vn_def, {'ROI_FullName_EN','Merged_ROI_EN'});
defRoiAbbrCol = pick_col_robust(vn_def, {'ROI_Abbrev','Merged_ROI_Abbrev'});
defRoiCnCol   = pick_col_robust(vn_def, {'ROI_CN','Merged_ROI_CN'});
if isempty(defRoiIdCol) || isempty(defRoiEnCol) || isempty(defRoiAbbrCol) || isempty(defRoiCnCol)
    error('ROI_Definition sheet must contain ROI_ID, ROI_FullName_EN, ROI_Abbrev and ROI_CN.');
end

hcpexIds = double(Tmap.(idCol));
roiIds   = double(Tmap.(roiIdCol));
valid = isfinite(hcpexIds) & isfinite(roiIds) & hcpexIds > 0 & roiIds > 0;
Tmap = Tmap(valid, :);
hcpexIds = double(Tmap.(idCol));
roiIds   = double(Tmap.(roiIdCol));

fprintf('\n[CHECK] Validating Parcel_Map...\n');

[uniqueIds, ~, idIC] = unique(hcpexIds);
idCounts = accumarray(idIC, 1);
dupIds = uniqueIds(idCounts > 1);
if ~isempty(dupIds)
    Tdup = Tmap(ismember(hcpexIds, dupIds), :);
    writetable(Tdup, fullfile(CFG.OutDir, 'ERROR_duplicate_HCPex_ID_assignments.csv'));
    error('Duplicated HCPex_ID assignments found. See ERROR_duplicate_HCPex_ID_assignments.csv');
end

defRoiIds = double(Tdef.(defRoiIdCol));
missingDef = setdiff(unique(roiIds), defRoiIds);
if ~isempty(missingDef)
    error('Some Merged_ROI_ID values in Parcel_Map are missing from ROI_Definition: %s', mat2str(missingDef(:)'));
end

if numel(unique(defRoiIds)) ~= numel(defRoiIds)
    error('ROI_Definition contains duplicated ROI_ID values.');
end

fprintf('[CHECK] Parcel_Map rows             : %d\n', height(Tmap));
fprintf('[CHECK] Unique HCPex IDs assigned   : %d\n', numel(uniqueIds));
fprintf('[CHECK] Merged ROI count in def     : %d\n', numel(unique(defRoiIds)));

fprintf('\n[BUILD] Reading HCPex atlas...\n');
V = spm_vol(CFG.HCPexAtlasNii);
Y = round(spm_read_vols(V));

atlasIds = unique(Y(:));
atlasIds = atlasIds(isfinite(atlasIds) & atlasIds > 0);
missingInAtlas = setdiff(uniqueIds, atlasIds);
if ~isempty(missingInAtlas)
    Tmiss = table(missingInAtlas(:), 'VariableNames', {'HCPex_ID_missing_in_atlas'});
    writetable(Tmiss, fullfile(CFG.OutDir, 'ERROR_HCPex_ID_missing_in_atlas.csv'));
    error('Some HCPex_ID values in Parcel_Map are not present in atlas. See ERROR_HCPex_ID_missing_in_atlas.csv');
end

outY = zeros(size(Y), 'double');
for i = 1:height(Tmap)
    sid = double(Tmap.(idCol)(i));
    mid = double(Tmap.(roiIdCol)(i));
    outY(Y == sid) = mid;
end

Vo = V;
Vo.fname = CFG.OutNii;
try
    Vo.dt = [spm_type('uint16') 0];
    Vo.pinfo = [1; 0; 0];
catch
    Vo.pinfo = [1; 0; 0];
end
spm_write_vol(Vo, outY);
fprintf('[BUILD] Wrote merged atlas: %s\n', CFG.OutNii);

fprintf('\n[QC] Writing lookup and voxel-count tables...\n');

Tlookup = table();
Tlookup.HCPex_ID          = double(Tmap.(idCol));
Tlookup.HCPex_LR          = cellstr(string(Tmap.(lrCol)));
Tlookup.HCPex_Parcel      = cellstr(string(Tmap.(parcelCol)));
Tlookup.Merged_ROI_ID     = double(Tmap.(roiIdCol));
Tlookup.Merged_ROI_Abbrev = cellstr(string(Tmap.(roiAbbrCol)));
Tlookup.Merged_ROI_EN     = cellstr(string(Tmap.(roiEnCol)));
Tlookup.Merged_ROI_CN     = cellstr(string(Tmap.(roiCnCol)));
Tlookup = sortrows(Tlookup, {'Merged_ROI_ID','HCPex_ID'});
writetable(Tlookup, CFG.OutLookupCsv);

defRoiIds = double(Tdef.(defRoiIdCol));
[defRoiIdsSorted, ord] = sort(defRoiIds, 'ascend');
TdefSorted = Tdef(ord, :);

Rows = {};
roi_ids = defRoiIdsSorted(:);
roi_labels = cell(numel(roi_ids), 1);
roi_abbrev = cell(numel(roi_ids), 1);
roi_cn = cell(numel(roi_ids), 1);
roi_voxels = zeros(numel(roi_ids), 1);

for r = 1:numel(roi_ids)
    rid = roi_ids(r);
    idxDef = find(double(TdefSorted.(defRoiIdCol)) == rid, 1, 'first');
    roiLabel = char(string(TdefSorted.(defRoiEnCol)(idxDef)));
    roiAbbr  = char(string(TdefSorted.(defRoiAbbrCol)(idxDef)));
    roiCN    = char(string(TdefSorted.(defRoiCnCol)(idxDef)));

    idxMap = find(double(Tlookup.Merged_ROI_ID) == rid);
    hids = Tlookup.HCPex_ID(idxMap);
    hparcels = unique(string(Tlookup.HCPex_Parcel(idxMap)), 'stable');
    nvox = sum(outY(:) == rid);

    roi_labels{r} = roiLabel;
    roi_abbrev{r} = roiAbbr;
    roi_cn{r} = roiCN;
    roi_voxels(r) = nvox;

    Rows(end+1,:) = {rid, roiAbbr, roiLabel, roiCN, numel(hparcels), numel(hids), nvox, ...
        strjoin(string(hids(:)'), ';'), strjoin(hparcels(:)', '; ')}; %#ok<AGROW>
end

Tvox = cell2table(Rows, 'VariableNames', ...
    {'ROI_ID','ROI_Abbrev','ROI_EN','ROI_CN','N_Parcels_Unilateral','N_HCPex_IDs_Bilateral','NVox_AtlasSpace','HCPex_IDs','HCPex_Parcels'});
writetable(Tvox, CFG.OutVoxelCsv);

if any(Tvox.NVox_AtlasSpace <= 0)
    warning('Some merged ROI has zero voxels. Please inspect roi_voxelcount.csv before continuing.');
end

save(CFG.OutCacheMAT, 'roi_ids', 'roi_labels', 'roi_abbrev', 'roi_cn', 'roi_voxels', 'Tlookup', 'Tvox');
fprintf('[QC] Wrote lookup csv     : %s\n', CFG.OutLookupCsv);
fprintf('[QC] Wrote voxelcount csv : %s\n', CFG.OutVoxelCsv);
fprintf('[QC] Wrote cache MAT      : %s\n', CFG.OutCacheMAT);

fprintf('\n[DONE] HCPex Merged17 atlas generation finished.\n');
fprintf('Please inspect NVox_AtlasSpace in:\n  %s\n', CFG.OutVoxelCsv);

end

function mkdir_if_needed(d)
if ~exist(d, 'dir')
    mkdir(d);
end
end

function T = readtable_preserve(fp, varargin)
try
    T = readtable(fp, 'VariableNamingRule', 'preserve', varargin{:});
catch
    T = readtable(fp, 'PreserveVariableNames', true, varargin{:});
end
end

function T = normalize_table_columns(T)
for i = 1:numel(T.Properties.VariableNames)
    col = T.Properties.VariableNames{i};
    if iscell(T.(col))
        T.(col) = cellfun(@normalize_cell_text, T.(col), 'UniformOutput', false);
    end
end
end

function y = normalize_cell_text(x)
if ischar(x) || isstring(x)
    y = char(strtrim(string(x)));
else
    y = x;
end
end

function col = pick_col_robust(vn, candidates)
col = '';
vn0 = standardize_names(vn);
for i = 1:numel(candidates)
    c0 = standardize_name_one(candidates{i});
    idx = find(strcmp(vn0, c0), 1, 'first');
    if ~isempty(idx)
        col = vn{idx};
        return;
    end
end
end

function out = standardize_names(x)
if iscell(x)
    out = cellfun(@standardize_name_one, x, 'UniformOutput', false);
else
    out = standardize_name_one(x);
end
end

function s = standardize_name_one(s)
s = char(string(s));
s = strtrim(s);
s = strrep(s, char(160), ' ');
s = regexprep(s, '\s+', '');
s = lower(s);
end
