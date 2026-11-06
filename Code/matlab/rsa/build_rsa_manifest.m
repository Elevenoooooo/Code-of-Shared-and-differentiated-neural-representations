function manifest = build_rsa_manifest(cfg)
%BUILD_RSA_MANIFEST Associate LSS beta files with stimulus identities/ratings.
% Run run_build_conditions and run_lss_first_level first. StimMeta_RT
% provides BetaIndex, OriginalTrial, ImageName and KeyRaw for all 99 trials.
required={'conditions_root','lss_root','derivatives_root','lss_run_manifest'};
assert(all(isfield(cfg,required)),'Missing RSA manifest-builder configuration.');
files=dir(fullfile(cfg.conditions_root,'**','StimMeta_RT_Group*_*.csv'));
assert(~isempty(files),'No StimMeta_RT sidecars. Run run_build_conditions first.');
out=fullfile(cfg.derivatives_root,'rsa','trial_tables'); if ~isfolder(out), mkdir(out); end
rows=cell(0,6); missing=cell(0,2);
taskMap=["Safety","Vitality","Beauty","Imageability"; ...
    "Vitality","Safety","Imageability","Beauty"; ...
    "Beauty","Imageability","Safety","Vitality"; ...
    "Imageability","Beauty","Vitality","Safety"];
for i=1:numel(files)
    token=regexp(files(i).name,'^StimMeta_RT_Group([1-4])_(\d{3})_Run([1-4])\.csv$','tokens','once');
    if isempty(token), continue; end
    g=str2double(token{1}); subject=token{2}; r=str2double(token{3});
    group="Group"+g; run="Run"+r; task=taskMap(g,r);
    betaDir=fullfile(cfg.lss_root,group,run,"UP_"+subject,"LSS_StimBeta_"+run);
    if ~isfolder(betaDir)
        missing(end+1,:)={string(files(i).name),string(betaDir)}; %#ok<AGROW>
        continue
    end
    source=fullfile(files(i).folder,files(i).name);
    T=readtable(source,'Delimiter',',','ReadVariableNames',true,'VariableNamingRule','preserve','TextType','string');
    assert(all(ismember({'BetaIndex','OriginalTrial','ImageName','KeyRaw'},T.Properties.VariableNames)), ...
        'StimMeta sidecar lacks BetaIndex,OriginalTrial,ImageName,KeyRaw.');
    assert(height(T)==99,'Each stimulus metadata file must contain 99 retained trials.');
    assert(isequal(double(T.BetaIndex),(1:height(T))'),'Noncontiguous stimulus metadata BetaIndex.');
    BetaFile=compose('beta_stim_trial_%03d.nii',double(T.BetaIndex));
    assert(all(isfile(fullfile(betaDir,BetaFile))),'A retained trial has no corresponding LSS beta.');
    count=dir(fullfile(betaDir,'beta_stim_trial_*.nii'));
    assert(numel(count)==height(T),'Beta count differs from stimulus metadata.');
    ImageName=rsa_image_ids(T.ImageName)+".png";
    Rating=parse_rating(T.KeyRaw); OriginalTrial=T.OriginalTrial; BetaIndex=T.BetaIndex;
    SourceStimMeta=repmat(string(source),height(T),1);
    trials=table(BetaFile,ImageName,Rating,BetaIndex,OriginalTrial,SourceStimMeta);
    trialFile=fullfile(out,sprintf('%s_%s_%s_%s_trials.csv',group,subject,run,task));
    writetable(trials,trialFile);
    rows(end+1,:)={group,string(subject),run,task,string(betaDir),string(trialFile)}; %#ok<AGROW>
end
manifest=cell2table(rows,'VariableNames',{'Group','Subject','Run','Task','BetaDir','TrialTableFile'});
assert(~isempty(manifest),'No LSS runs have matching stimulus metadata.');
parent=fileparts(cfg.lss_run_manifest); if ~isfolder(parent), mkdir(parent); end
writetable(manifest,cfg.lss_run_manifest);
writetable(cell2table(missing,'VariableNames',{'StimMetaFile','MissingBetaDir'}),fullfile(out,'missing_lss_runs.csv'));
end

function x=parse_rating(v)
s=lower(strtrim(string(v))); s=erase(s,'num_'); s=erase(s,'key_'); s=erase(s,'button_');
s=regexprep(s,'^[''" ]+|[''" ]+$',''); x=str2double(s);
for i=find(~isfinite(x))'
    token=regexp(char(s(i)),'[-+]?\d+(\.\d+)?','match','once');
    if ~isempty(token), x(i)=str2double(token); end
end
x(ismember(s,["yes","y","true","positive"]))=2;
x(ismember(s,["no","n","false","negative"]))=1;
end
