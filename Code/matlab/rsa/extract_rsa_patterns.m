function runManifest = extract_rsa_patterns(cfg)
%EXTRACT_RSA_PATTERNS Extract ROI voxel patterns from single-trial LSS NIfTIs.
% Required cfg.lss_run_manifest CSV: Group,Subject,Run,Task,BetaDir,TrialTableFile.
% Trial tables explicitly pair each BetaFile with ImageName and Rating.
% This function samples merged atlas labels in each run's native beta grid
% with nearest-neighbour interpolation and retains each ROI's voxel patterns.
assert(isfield(cfg,'spm_dir')&&isfolder(cfg.spm_dir),'Set cfg.spm_dir to SPM12.');
addpath(cfg.spm_dir); spm('Defaults','fMRI');
assert(isfile(cfg.hcpex_atlas)&&isfile(cfg.hcpex_labels),'Missing merged HCPex atlas or labels.');
M=readtable(cfg.lss_run_manifest,'Delimiter',',','ReadVariableNames',true,'VariableNamingRule','preserve','TextType','string');
assert(all(ismember({'Group','Subject','Run','Task','BetaDir','TrialTableFile'},M.Properties.VariableNames)), ...
    'LSS manifest schema mismatch.');
L=readtable(cfg.hcpex_labels,'Delimiter',',','ReadVariableNames',true,'VariableNamingRule','preserve','TextType','string');
if all(ismember({'Merged_ROI_ID','Merged_ROI_EN'},L.Properties.VariableNames))
    [roi_ids,ia]=unique(double(L.Merged_ROI_ID),'stable'); roi_labels=cellstr(L.Merged_ROI_EN(ia));
else
    assert(all(ismember({'ROI_ID','Name'},L.Properties.VariableNames)), ...
        'Labels require Merged_ROI_ID,Merged_ROI_EN or packaged ROI_ID,Name.');
    [roi_ids,ia]=unique(double(L.ROI_ID),'stable'); roi_labels=cellstr(L.Name(ia));
end
assert(numel(roi_ids)==17,'Expected the study merged bilateral HCPex17 atlas.');
VAtlas=spm_vol(cfg.hcpex_atlas); assert(numel(VAtlas)==1,'Atlas must be a single 3D label map.');
root=fullfile(cfg.derivatives_root,'rsa','run_patterns'); if ~isfolder(root), mkdir(root); end
rows=cell(0,6);
for r=1:height(M)
    group=char(M.Group(r)); subject=sprintf('%03d',str2double(string(M.Subject(r))));
    run=char(M.Run(r)); task=char(M.Task(r));
    assert(ismember(task,{'Safety','Vitality','Beauty','Imageability'}),'Unknown task.');
    betaDir=resolve_path(M.BetaDir(r),cfg.derivatives_root,false);
    trialFile=resolve_path(M.TrialTableFile(r),cfg.derivatives_root,true);
    T=readtable(trialFile,'Delimiter',',','ReadVariableNames',true,'VariableNamingRule','preserve','TextType','string');
    assert(all(ismember({'BetaFile','ImageName','Rating'},T.Properties.VariableNames)), ...
        'Trial table needs BetaFile,ImageName,Rating.');
    assert(height(T)==99,'Each trial table must contain 99 retained stimulus trials.');
    imageIDs=rsa_image_ids(T.ImageName);
    assert(numel(unique(imageIDs))==height(T),'Repeated images in the run trial table.');
    betas=cell(height(T),1);
    for t=1:height(T), betas{t}=resolve_path(T.BetaFile(t),betaDir,true); end
    assert(numel(unique(string(betas)))==height(T),'Duplicate BetaFile entries.');
    Vref=spm_vol(betas{1}); assert(numel(Vref)==1,'Each beta must be a 3D NIfTI.');
    [i,j,k]=ndgrid(1:Vref.dim(1),1:Vref.dim(2),1:Vref.dim(3));
    locations=VAtlas.mat\(Vref.mat*[i(:)';j(:)';k(:)';ones(1,numel(i))]);
    atlas=spm_sample_vol(VAtlas,locations(1,:),locations(2,:),locations(3,:),0);
    atlas=round(atlas(:));
    indices=cell(numel(roi_ids),1); PatternCell=cell(size(indices));
    for a=1:numel(roi_ids)
        indices{a}=find(atlas==roi_ids(a));
        PatternCell{a}=nan(height(T),numel(indices{a}));
    end
    for t=1:height(T)
        V=spm_vol(betas{t});
        assert(numel(V)==1&&isequal(V.dim,Vref.dim)&&max(abs(V.mat-Vref.mat),[],'all')<1e-5, ...
            'Beta geometries differ within a run.');
        data=spm_read_vols(V);
        for a=1:numel(roi_ids), PatternCell{a}(t,:)=reshape(data(indices{a}),1,[]); end
    end
    pathRaw=imageIDs+".png"; ratingVec=str2double(string(T.Rating));
    GroupName=string(group); Subject=string(subject); RunName=string(run); TaskName=string(task);
    SourceTrialTable=string(trialFile); BetaFiles=string(betas);
    nBeta=height(T); SourceAtlas=string(cfg.hcpex_atlas);
    destination=fullfile(root,sprintf('UP_%s_%s_%s_%s_HCPex17_RunPattern.mat',group,subject,run,task));
    if ismember('RT',T.Properties.VariableNames)
        rtVec=str2double(string(T.RT));
    else
        rtVec=nan(nBeta,1);
    end
    save(destination,'PatternCell','pathRaw','ratingVec','rtVec','roi_ids','roi_labels', ...
        'GroupName','Subject','RunName','TaskName','nBeta','SourceTrialTable','BetaFiles','SourceAtlas','-v7.3');
    rows(end+1,:)={GroupName,Subject,RunName,TaskName,string(destination),string(trialFile)}; %#ok<AGROW>
end
runManifest=cell2table(rows,'VariableNames',{'Group','Subject','Run','Task','RunPatternFile','RatingTableFile'});
parent=fileparts(cfg.rsa_run_manifest); if ~isfolder(parent), mkdir(parent); end
writetable(runManifest,cfg.rsa_run_manifest);
end

function path=resolve_path(value,root,isFile)
path=char(value);
if isFile
    if ~isfile(path), path=fullfile(root,path); end
    assert(isfile(path),'Missing file: %s',path);
else
    if ~isfolder(path), path=fullfile(root,path); end
    assert(isfolder(path),'Missing folder: %s',path);
end
end
