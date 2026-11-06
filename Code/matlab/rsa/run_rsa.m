function results = run_rsa(cfg)
%RUN_RSA Figure 5 environmental-feature and perceptual-rating RSA.
% Required cfg: derivatives_root,output_root,
% urban_rdm_dir,rsa_run_manifest. See docs/rsa.md for manifest schemas.
% Each neural and rating RDM represents 99 retained stimulus trials.
cfg = defaults(cfg);
out = fullfile(cfg.output_root, 'rsa'); if ~isfolder(out), mkdir(out); end
M = readtable(cfg.rsa_run_manifest, 'Delimiter',',','ReadVariableNames',true,'VariableNamingRule','preserve', 'TextType','string');
required = {'Group','Subject','Run','Task','RunPatternFile'};
assert(all(ismember(required,M.Properties.VariableNames)), 'Invalid RSA manifest columns: %s',strjoin(M.Properties.VariableNames,','));
assert(height(M)>0, 'RSA manifest is empty.');
tasks = ["Safety","Vitality","Beauty","Imageability"];
assert(all(ismember(M.Task,tasks)), 'Manifest has unknown Task values.');
M.Subject = compose('%03d', str2double(string(M.Subject)));
assert(all(~contains(M.Subject,'NaN')), 'Subject must have a numeric study ID.');
M.SubjectKey = M.Group + "_" + M.Subject;
key = M.SubjectKey + "_" + M.Run + "_" + M.Task;
assert(numel(unique(key)) == height(M), 'Duplicate subject/run/task rows in manifest.');
keep = true(height(M),1);
if cfg.rsa_complete_cases
    subjects = unique(M.SubjectKey);
    for i=1:numel(subjects)
        rows = M.SubjectKey==subjects(i);
        keep(rows) = all(ismember(tasks,unique(M.Task(rows))));
    end
end
M.IncludedCompleteCase = keep;
writetable(M, fullfile(out,'run_manifest_inclusion.csv'));
M = M(keep,:); assert(height(M)>0, 'No complete-case subjects in manifest.');
[predNames, predRDMs] = load_features(cfg);
pair = tril(true(400),-1);
Xglobal = zeros(sum(pair(:)),10);
for p=1:10, current=predRDMs(:,:,p); Xglobal(:,p)=current(pair); end
QCcorr = array2table(corr(Xglobal),'VariableNames',cellstr(predNames),'RowNames',cellstr(predNames));
writetable(QCcorr,fullfile(out,'predictor_rdm_pearson.csv'),'WriteRowNames',true);
writetable(table(predNames, vif(zstandard(Xglobal)), 'VariableNames',{'Predictor','VIF'}), ...
    fullfile(out,'predictor_rdm_vif.csv'));
betas=cell(0,11); ratings=cell(0,10); qc=cell(0,7); fits=cell(0,7);
for f=1:height(M)
    path = resolve_file(M.RunPatternFile(f),cfg.derivatives_root);
    S = load(path);
    assert(isfield(S,'PatternCell') && isfield(S,'roi_labels'), 'Missing PatternCell/roi_labels: %s',path);
    ids = run_ids(S,cfg);
    assert(numel(unique(ids))==numel(ids), 'Repeated images in run: %s',path);
    nTrial = numel(ids);
    if cfg.rsa_require_expected_trials
        assert(nTrial==cfg.rsa_expected_trials,'Expected %d retained trials in each run.',cfg.rsa_expected_trials);
    end
    rating = run_ratings(S,M(f,:),ids,cfg);
    labels = string(S.roi_labels(:));
    lowerLabels = lower(labels);
    keepROI = ~contains(lowerLabels,'amyg') & ~contains(lowerLabels,'entorhinal') & ~strcmpi(labels,'EC');
    assert(sum(keepROI)==15, 'Expected 15 retained ROIs after Amygdala/Entorhinal exclusion: %s',path);
    [found, imageIndex] = ismember(ids,"image"+string((1:400)'));
    assert(all(found), 'Unknown canonical ID.');
    selected = select_three(labels);
    assert(numel(selected)==3, 'Could not identify exactly VVC/PPC/vmPFC.');
    for ri=find(keepROI)'
        X = S.PatternCell{ri};
        assert(size(X,1)==nTrial,'Trial count differs between pathRaw and PatternCell.');
        try
            [D, trialKeep, nVox] = rsa_neural_rdm(X,cfg.rsa_min_trials,cfg.rsa_min_voxels);
        catch ME
            qc(end+1,:)={M.SubjectKey(f),M.Task(f),labels(ri),nTrial,0,0,string(ME.message)}; %#ok<AGROW>
            continue
        end
        if cfg.rsa_require_expected_trials
            assert(size(D,1)==cfg.rsa_expected_trials,'Expected %d finite trial patterns in each neural RDM.',cfg.rsa_expected_trials);
        end
        qc(end+1,:)={M.SubjectKey(f),M.Task(f),labels(ri),nTrial,sum(trialKeep),nVox,"ok"}; %#ok<AGROW>
        % Every response and predictor uses the same lower-triangle order.
        mask=tril(true(size(D)),-1); y=D(mask);
        R=abs(rating(trialKeep)-rating(trialKeep)'); x=R(mask);
        ok=isfinite(x)&isfinite(y);
        if sum(ok)>=cfg.rsa_min_pairs
            rho=corr(x(ok),y(ok),'Type','Spearman','Rows','complete');
            fisherZ=atanh(max(min(rho,0.999999),-0.999999));
            ratings(end+1,:)={M.SubjectKey(f),M.Group(f),M.Subject(f),M.Run(f),M.Task(f),ri,labels(ri),sum(ok),rho,fisherZ}; %#ok<AGROW>
        end
        if ~ismember(ri,selected), continue; end
        index=imageIndex(trialKeep); Xpred=zeros(numel(y),10);
        for p=1:10, R=predRDMs(index,index,p); Xpred(:,p)=R(mask); end
        valid=isfinite(y)&all(isfinite(Xpred),2);
        if sum(valid)<cfg.rsa_min_pairs, continue; end
        Xpred=Xpred(valid,:); y=y(valid);
        constant=std(Xpred,0,1)<=eps;
        if all(constant) || std(y)<=eps, continue; end
        [b, fit, vf, uniqueR2]=fit_model(zstandard(Xpred(:,~constant)),zstandard(y));
        positions=find(~constant);
        for p=1:numel(positions)
            betas(end+1,:)={M.SubjectKey(f),M.Group(f),M.Subject(f),M.Run(f),M.Task(f),ri,labels(ri),predNames(positions(p)),b(p),vf(p),uniqueR2(p)}; %#ok<AGROW>
        end
        fits(end+1,:)={M.SubjectKey(f),M.Task(f),labels(ri),sum(valid),fit.R2,fit.AdjR2,strjoin(predNames(constant),';')}; %#ok<AGROW>
    end
end
Tbeta=cell2table(betas,'VariableNames',{'SubjectKey','GroupName','Subject','Run','Task','ROI_ID','ROI_Label','Predictor','BetaStd','VIF','UniqueR2'});
Trating=cell2table(ratings,'VariableNames',{'SubjectKey','GroupName','Subject','Run','Task','ROI_ID','ROI_Label','N_Pairs','SpearmanRho','FisherZ'});
Tqc=cell2table(qc,'VariableNames',{'SubjectKey','Task','ROI_Label','N_InputTrials','N_TrialUsed','N_VoxelKept','Status'});
Tfit=cell2table(fits,'VariableNames',{'SubjectKey','Task','ROI_Label','N_Pairs','R2','AdjR2','DroppedConstantPredictors'});
assert(~isempty(Tbeta)&&~isempty(Trating),'No analyzable RSA output; inspect manifest and voxel/trial QC.');
TbetaStat=group_inference(Tbeta,'BetaStd',cfg.rsa_min_group_n,false,cfg);
TratingStat=group_inference(Trating,'FisherZ',2,true,cfg);
writetable(Tbeta,fullfile(out,'subject_run_feature_beta.csv'));
writetable(Trating,fullfile(out,'subject_run_rating_rsa.csv'));
writetable(TbetaStat,fullfile(out,'group_feature_beta.csv'));
writetable(TratingStat,fullfile(out,'group_rating_rsa.csv'));
writetable(Tqc,fullfile(out,'neural_rdm_qc.csv'));
writetable(Tfit,fullfile(out,'model_fit.csv'));
results=struct('feature_beta',Tbeta,'rating_rsa',Trating,'feature_group',TbetaStat, ...
    'rating_group',TratingStat,'neural_qc',Tqc,'model_fit',Tfit);
save(fullfile(out,'rsa_results.mat'),'results','cfg','-v7.3');
end

function cfg=defaults(cfg)
fields={'rsa_expected_trials','rsa_require_expected_trials','rsa_min_trials','rsa_min_voxels','rsa_min_pairs','rsa_min_group_n','rsa_complete_cases'};
values={99,true,20,10,1000,8,true};
for i=1:numel(fields), if ~isfield(cfg,fields{i}), cfg.(fields{i})=values{i}; end, end
if ~isfield(cfg,'rsa_feature_fdr_mode'), cfg.rsa_feature_fdr_mode='SelectedROI_Global'; end
end

function [names,rdms]=load_features(cfg)
names=["Colorfulness";"Visual Complexity";"Spatial Layering";"Green View Index";"Sky Visibility";"Enclosure";"Crowd Attraction";"Vehicle Interference";"Street Facility Convenience";"Urbanization Intensity Index"];
files=["colorfulness";"visual_complexity";"spatial_layering";"green_view_index";"sky_visibility";"enclosure";"crowd_attraction";"vehicle_interference";"street_facility_convenience";"urbanization_intensity_index"];
order="image"+string((1:400)'); rdms=zeros(400,400,10);
for p=1:10
    fp=fullfile(cfg.urban_rdm_dir,"RDM_"+files(p)+".csv");
    T=readtable(fp,'Delimiter',',','ReadVariableNames',true,'VariableNamingRule','preserve','TextType','string');
    row=rsa_image_ids(T{:,1}); col=rsa_image_ids(T.Properties.VariableNames(2:end));
    assert(numel(row)==400&&numel(col)==400&&numel(unique(row))==400&&numel(unique(col))==400,'RDM must have 400 unique labels.');
    [okR,ir]=ismember(order,row); [okC,ic]=ismember(order,col);
    assert(all(okR)&all(okC),'Predictor RDM is missing stimulus labels.');
    R=double(T{:,2:end}); R=R(ir,ic);
    assert(all(isfinite(R(:)))&&all(R(:)>=-1e-8),'Nonfinite/negative objective-feature distance.');
    assert(max(abs(R-R'),[],'all')<=1e-8&&max(abs(diag(R)))<=1e-8,'RDM symmetry/diagonal mismatch.');
    rdms(:,:,p)=R;
end
end

function ids=run_ids(S,cfg)
if isfield(S,'pathRaw'), ids=rsa_image_ids(S.pathRaw); return; end
if isfield(S,'ImageNames'), ids=rsa_image_ids(S.ImageNames); return; end
assert(isfield(S,'canonIdx')&&isfield(cfg,'image_id_csv'), ...
    'Run patterns need pathRaw/ImageNames or canonIdx with cfg.image_id_csv.');
T=readtable(cfg.image_id_csv,'Delimiter',',','ReadVariableNames',true,'VariableNamingRule','preserve','TextType','string');
assert(all(ismember({'File_Name','Img_ID'},T.Properties.VariableNames)),'Image manifest requires File_Name,Img_ID.');
experimentIndex=str2double(string(T.Img_ID));
assert(all(isfinite(experimentIndex))&&numel(unique(experimentIndex))==height(T),'Experiment indices must be finite and unique.');
canonical=rsa_image_ids(T.File_Name);
assert(numel(canonical)==400&&numel(unique(canonical))==400,'Canonical file must specify all 400 unique images.');
% canonIdx values are experiment indices joined to Img_ID, not filename numbers.
idx=double(S.canonIdx(:)); [found,positions]=ismember(idx,experimentIndex);
assert(all(found),'A canonIdx experiment index is absent from the image manifest.');
ids=canonical(positions);
end

function values=run_ratings(S,row,ids,cfg)
if isfield(S,'ratingVec')&&numel(S.ratingVec)==numel(ids)
    values=double(S.ratingVec(:)); return
end
assert(ismember('RatingTableFile',row.Properties.VariableNames), ...
    'Missing ratingVec. Provide manifest RatingTableFile with explicit ImageName,Rating in beta order.');
T=readtable(resolve_file(row.RatingTableFile,cfg.derivatives_root),'Delimiter',',','ReadVariableNames',true,'VariableNamingRule','preserve','TextType','string');
assert(all(ismember({'ImageName','Rating'},T.Properties.VariableNames)),'Rating table needs ImageName,Rating.');
keys=rsa_image_ids(T.ImageName); assert(numel(unique(keys))==numel(keys),'Duplicate rating labels.');
[found,pos]=ismember(ids,keys); assert(all(found),'Rating table misses beta images.');
values=str2double(string(T.Rating(pos)));
end

function ids=select_three(labels)
low=lower(labels);
a=find(strcmpi(labels,'VVC')|contains(low,'ventral visual cortex'));
b=find(strcmpi(labels,'PPC')|strcmpi(labels,'SPL/IPL')|contains(low,'superior and inferior parietal lobules'));
c=find(strcmpi(labels,'vmPFC')|contains(low,'ventromedial prefrontal cortex'));
assert(numel(a)==1&&numel(b)==1&&numel(c)==1,'Ambiguous VVC/PPC/vmPFC ROI labels.');
ids=[a;b;c];
end

function [b,fit,vf,u]=fit_model(X,y)
D=[ones(numel(y),1),X]; beta=D\y; b=beta(2:end);
sst=sum((y-mean(y)).^2); sse=sum((y-D*beta).^2);
fit.R2=1-sse/max(sst,eps); fit.AdjR2=1-(1-fit.R2)*(numel(y)-1)/max(numel(y)-size(X,2)-1,1);
vf=vif(X); u=zeros(size(b));
for p=1:size(X,2)
    red=D; red(:,p+1)=[]; r2red=1-sum((y-red*(red\y)).^2)/max(sst,eps);
    u(p)=fit.R2-r2red;
end
end

function result=vif(X)
result=ones(size(X,2),1);
for p=1:size(X,2)
    other=X; other(:,p)=[];
    if isempty(other), continue; end
    D=[ones(size(X,1),1),other]; y=X(:,p); sst=sum((y-mean(y)).^2);
    r2=1-sum((y-D*(D\y)).^2)/max(sst,eps);
    result(p)=1/max(1-r2,eps);
end
end

function X=zstandard(X)
sd=std(X,0,1); sd(~isfinite(sd)|sd<=0)=1; X=(X-mean(X,1))./sd;
end

function T=group_inference(source,value,minN,isRating,cfg)
if isRating, source.Predictor=repmat("rating",height(source),1); end
keys=unique(source(:,{'Task','ROI_ID','ROI_Label','Predictor'}),'rows','stable'); rows={};
for k=1:height(keys)
    take=source.Task==keys.Task(k)&source.ROI_ID==keys.ROI_ID(k)&source.Predictor==keys.Predictor(k);
    subset=source(take,:); [g,~]=findgroups(subset.SubjectKey);
    x=splitapply(@(v)mean(v,'omitnan'),subset.(value),g); x=x(isfinite(x)); n=numel(x);
    sd=std(x); mu=mean(x); p=NaN; t=NaN; df=NaN;
    if n>=minN&&isfinite(sd)&&sd>eps
        [~,p,~,st]=ttest(x,0,'Tail','both'); t=st.tstat; df=st.df;
    end
    meanRho=NaN;
    if isRating
        rho=splitapply(@(v)mean(v,'omitnan'),subset.SpearmanRho,g); meanRho=mean(rho,'omitnan');
    end
    rows(end+1,:)={keys.Task(k),keys.ROI_ID(k),keys.ROI_Label(k),keys.Predictor(k),n,mu,sd,sd/sqrt(max(n,1)),t,df,p,meanRho}; %#ok<AGROW>
end
T=cell2table(rows,'VariableNames',{'Task','ROI_ID','ROI_Label','Predictor','N_Subjects','Mean','SD','SEM','T','DF','P_uncorrected','MeanSpearmanRho'});
T.Q_FDR_Global=bh(T.P_uncorrected); T.Q_FDR_SelectedROI_Global=T.Q_FDR_Global;
T.Q_FDR_WithinROI=nan(height(T),1); T.Q_FDR_WithinTask=nan(height(T),1);
T.Q_FDR_WithinTaskROI=nan(height(T),1); T.FDR_FamilySize=zeros(height(T),1);
rois=unique(T.ROI_ID);
for i=1:numel(rois)
    idx=T.ROI_ID==rois(i); T.Q_FDR_WithinROI(idx)=bh(T.P_uncorrected(idx));
    if isRating, T.FDR_FamilySize(idx)=sum(isfinite(T.P_uncorrected(idx))); end
end
tasks=unique(T.Task);
for i=1:numel(tasks)
    idx=T.Task==tasks(i); T.Q_FDR_WithinTask(idx)=bh(T.P_uncorrected(idx));
    for j=1:numel(rois)
        sub=idx&T.ROI_ID==rois(j); T.Q_FDR_WithinTaskROI(sub)=bh(T.P_uncorrected(sub));
    end
end
if isRating
    T.Q_FDR=T.Q_FDR_WithinROI;
    T.FDR_Family="RatingWithinROI_AcrossTasks|ROI="+T.ROI_Label;
else
    mode=char(cfg.rsa_feature_fdr_mode);
    if strcmpi(mode,'Uncorrected')
        T.Q_FDR=T.P_uncorrected;
    else
        scope=['Q_FDR_' mode]; assert(ismember(scope,T.Properties.VariableNames),'Unknown rsa_feature_fdr_mode.');
        T.Q_FDR=T.(scope);
    end
    for k=1:height(T)
        switch lower(mode)
            case {'global','selectedroi_global'}, family=true(height(T),1);
            case 'withinroi', family=T.ROI_ID==T.ROI_ID(k);
            case 'withintask', family=T.Task==T.Task(k);
            case 'withintaskroi', family=T.Task==T.Task(k)&T.ROI_ID==T.ROI_ID(k);
            case 'uncorrected', family=false(height(T),1); family(k)=true;
        end
        T.FDR_FamilySize(k)=sum(isfinite(T.P_uncorrected(family)));
    end
    T.FDR_Family=repmat(string(mode),height(T),1);
end
T.Significant_FDR_0p05=T.Q_FDR<0.05;
end

function q=bh(p)
q=nan(size(p)); ok=isfinite(p); x=p(ok); [a,ord]=sort(x); n=numel(a);
if n==0, return; end
b=min(1,flipud(cummin(flipud(a.*n./(1:n)')))); result=nan(size(x)); result(ord)=b; q(ok)=result;
end

function p=resolve_file(value,root)
p=char(value); if ~isfile(p), p=fullfile(root,p); end
assert(isfile(p),'Missing input file: %s',p);
end
