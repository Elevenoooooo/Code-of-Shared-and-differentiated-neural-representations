function Statistics = fit_connectivity_task_lme(cfg,EdgeLong)
% FIT_CONNECTIVITY_TASK_LME Test task-related differences in 105 connections.
% Model: Z ~ Task + Group + (1|Subject), fitted by REML with effect coding.
% Omnibus task F tests use residual denominator degrees of freedom; BH-FDR
% correction is applied across the 105 connection tests at q<0.05.
% Accepts a supplied edge table or reads BSC_EdgeLong.csv.
if nargin < 1, cfg = study_config(); end
assert(exist('fitlme','file')==2,'Statistics and Machine Learning Toolbox is required.');
if nargin < 2
    EdgeLong = readtable(fullfile(cfg.connectivity_root,'BSC_EdgeLong.csv'), ...
        'TextType','string','VariableNamingRule','preserve');
end
required = {'Subject','Group','Task','Edge','Z','ROI1_ID','ROI1','ROI2_ID','ROI2'};
assert(all(ismember(required,EdgeLong.Properties.VariableNames)),'Missing edge table columns.');
opt = struct('FitMethod','REML','DFMethod','residual','MinimumSubjects',8);
if isfield(cfg,'connectivity_lme_options')
    fields = fieldnames(cfg.connectivity_lme_options);
    for i = 1:numel(fields)
        assert(isfield(opt,fields{i}),'Unknown connectivity_lme_options key: %s',fields{i});
        opt.(fields{i}) = cfg.connectivity_lme_options.(fields{i});
    end
end
edges = unique(EdgeLong.Edge); assert(numel(edges)==105,'Expected 105 edges for 15 ROIs.');
rows = cell(numel(edges),1);
taskNames = {'Safety','Vitality','Beauty','Imageability'};
for i = 1:numel(edges)
    T = EdgeLong(EdgeLong.Edge==edges(i),:);
    first = T(1,:);
    T.Subject = categorical(string(T.Subject));
    T.Group = categorical(string(T.Group));
    T.Task = categorical(string(T.Task),taskNames);
    T = T(isfinite(T.Z) & ~isundefined(T.Subject) & ~isundefined(T.Group) & ~isundefined(T.Task),:);
    subjects = categories(removecats(T.Subject)); include = false(height(T),1);
    for s = 1:numel(subjects)
        idx = T.Subject==subjects{s};
        assert(numel(unique(T.Task(idx)))==sum(idx),'Duplicate subject/task data for edge %d.',edges(i));
        if sum(idx)==4, include = include | idx; end
    end
    T = T(include,:); n = numel(unique(T.Subject));
    F = NaN; P = NaN; df1 = NaN; df2 = NaN; status = "Insufficient complete subjects";
    if n>=opt.MinimumSubjects
        assert(numel(unique(T.Group))==4,'All four Latin-square groups are required for edge %d.',edges(i));
        try
            lme = fitlme(T,'Z ~ Task + Group + (1|Subject)', ...
                'FitMethod',opt.FitMethod,'DummyVarCoding','effects');
            A = anova(lme,'DFMethod',opt.DFMethod);
            if isa(A,'dataset'), A = dataset2table(A); end
            row = strcmp(string(A.Term),'Task');
            assert(sum(row)==1,'Task effect missing from LME ANOVA.');
            F = A.FStat(row); P = A.pValue(row); df1 = A.DF1(row); df2 = A.DF2(row);
            status = "OK";
        catch ME
            status = "Failed: "+string(ME.message);
        end
    end
    rows{i} = table(edges(i),first.ROI1_ID,string(first.ROI1),first.ROI2_ID, ...
        string(first.ROI2),n,F,df1,df2,P,status, ...
        'VariableNames',{'Edge','ROI1_ID','ROI1','ROI2_ID','ROI2','NSubjects', ...
        'F','DF1','DF2','P','Status'});
end
Statistics = vertcat(rows{:});
Statistics.Q = bh_q(Statistics.P);
Statistics.Significant = Statistics.Q<0.05;
Statistics.Model = repmat("Z ~ Task + Group + (1|Subject)",height(Statistics),1);
Statistics.FitMethod = repmat(string(opt.FitMethod),height(Statistics),1);
Statistics.DFMethod = repmat(string(opt.DFMethod),height(Statistics),1);
if ~isfolder(cfg.connectivity_root), mkdir(cfg.connectivity_root); end
writetable(Statistics,fullfile(cfg.connectivity_root,'BSC_TaskMainEffect_LME.csv'));
assert(all(Statistics.Status=="OK"), ...
    'Some edge models were unavailable or failed. Inspect the saved Status column; do not interpret them as null results.');
end

function q = bh_q(p)
q = nan(size(p)); valid = isfinite(p) & p>=0 & p<=1;
[sorted,order] = sort(p(valid)); n = numel(sorted);
adjusted = min(1,flipud(cummin(flipud(sorted(:).*n./(1:n)'))));
back = zeros(n,1); back(order) = adjusted; q(valid) = back;
end
