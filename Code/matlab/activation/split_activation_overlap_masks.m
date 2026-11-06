function split_activation_overlap_masks(cfg)
% SPLIT_ACTIVATION_OVERLAP_MASKS Export masks for each task-overlap bit code.
% Reads the conjunction label map and writes 15 mutually exclusive masks and
% a legend recording the contributing tasks and voxel counts.

if nargin < 1, cfg = study_config(); end
required = {'spm_dir','activation_root'};
for k = 1:numel(required)
    assert(isfield(cfg,required{k}) && ~isempty(cfg.(required{k})), 'Missing cfg.%s', required{k});
end
addpath(cfg.spm_dir);

spm('Defaults','fMRI');

conjunctionDir = fullfile(cfg.activation_root,'FourTask_ActivationConjunction','Stim_gt_Fixation');
labelPath = fullfile(conjunctionDir,'Stim_gt_Fixation_FourTask_TaskOverlap_LabelMap_BitCode.nii');

outDir = fullfile(conjunctionDir,'BitCode_SplitMasks_ForMRIcroGL');

tasks = {'Safety','Vitality','Beauty','Imageability'};

if exist(outDir, 'dir') ~= 7
    mkdir(outDir);
end

V = spm_vol(labelPath);
Y = round(spm_read_vols(V));
Y(~isfinite(Y)) = 0;

rows = {};

for code = 1:15
    flags = false(1,4);
    for ti = 1:4
        flags(ti) = bitget(uint16(code), ti) == 1;
    end

    comboName = strjoin(tasks(flags), '_AND_');
    if isempty(comboName)
        comboName = 'None';
    end

    mask = uint8(Y == code);
    nVox = nnz(mask);

    outName = sprintf('BitCode_%02d_%s_mask.nii', code, comboName);
    outName = regexprep(outName, '[^\w\-\.\(\)]', '_');

    Vout = V;
    Vout.fname = fullfile(outDir, outName);
    Vout.dt = [2 0];
    spm_write_vol(Vout, mask);

    rows(end+1,:) = {code, comboName, sum(flags), nVox, Vout.fname}; %#ok<AGROW>

    fprintf('[WRITE] code=%02d | nVox=%d | %s\n', code, nVox, Vout.fname);
end

T = cell2table(rows, 'VariableNames', { ...
    'BitCode','TaskCombination','NTasks','VoxelCount','MaskPath'});

outCsv = fullfile(outDir, 'BitCode_Label_Legend.csv');
write_table_utf8bom(T, outCsv);
fprintf('[WRITE] %s\n', outCsv);

fprintf('\n[ALL DONE]\n');

end

function write_table_utf8bom(T, outPath)
fid = fopen(outPath, 'w', 'n', 'UTF-8');
fwrite(fid, char(65279), 'char');
fclose(fid);
writetable(T, outPath, 'FileType','text', 'Delimiter', ',', 'WriteMode','append');
end
