function T = make_rsa_trial_table(~, behaviorCSV, betaDir, behaviorOffset, outputCSV, ratingColumn)
%MAKE_RSA_TRIAL_TABLE Associate 99 LSS betas with behavior images and ratings.
% behaviorOffset is the number of initial formal trials omitted from LSS.
% The trial table preserves the selected behavior row and trial indices.
assert(nargin>=6&&~isempty(ratingColumn), ...
    'Supply the rating column and behaviorOffset.');
assert(isscalar(behaviorOffset)&&isfinite(behaviorOffset)&&behaviorOffset>=0&&mod(behaviorOffset,1)==0, ...
    'behaviorOffset must be a nonnegative integer.');
B=readtable(behaviorCSV,'Delimiter',',','ReadVariableNames',true,'VariableNamingRule','preserve','TextType','string');
assert(ismember(ratingColumn,B.Properties.VariableNames),'Missing requested rating column.');
imageColumns={'ImageName','Path_raw','Path'};
col='';
for i=1:numel(imageColumns), if ismember(imageColumns{i},B.Properties.VariableNames), col=imageColumns{i}; break; end, end
assert(~isempty(col),'No recognized image-path column in behavior CSV.');
% Select formal rows by standard image labels in the 400-stimulus set.
labels=basenames(B.(col)); formal=ismember(labels,"image"+string((1:400)'));
formalRows=find(formal); assert(numel(formalRows)==100,'Expected exactly 100 formal stimulus rows before offset.');
files=dir(fullfile(betaDir,'beta_stim_trial_*.nii'));
assert(numel(files)==99,'Expected 99 LSS beta_stim_trial_*.nii files.');
numbers=nan(numel(files),1);
for i=1:numel(files)
    token=regexp(files(i).name,'trial_(\d+)\.nii$','tokens','once'); numbers(i)=str2double(token{1});
end
[numbers,order]=sort(numbers); files=files(order);
assert(isequal(numbers,(1:numel(files))'),'LSS beta trial numbers must be contiguous, starting at 1.');
indices=(1:numel(files))'+behaviorOffset;
assert(max(indices)<=numel(formalRows),'Offset plus beta count exceeds available formal trials.');
selected=formalRows(indices);
ImageName=rsa_image_ids(B.(col)(selected))+".png";
Rating=parse_ratings(B.(ratingColumn)(selected));
BetaFile=string({files.name})'; BehaviorRow=selected; BehaviorTrial=indices;
BehaviorOffset=repmat(behaviorOffset,numel(files),1);
T=table(BetaFile,ImageName,Rating,BehaviorRow,BehaviorTrial,BehaviorOffset);
parent=fileparts(outputCSV); if ~isempty(parent)&&~isfolder(parent), mkdir(parent); end
writetable(T,outputCSV);
end

function result=basenames(values)
values=string(values(:)); result=strings(size(values));
for i=1:numel(values)
    [~,name,~]=fileparts(strrep(char(values(i)),'\','/')); result(i)=lower(string(name));
end
end

function x=parse_ratings(values)
% Convert numeric response values and key labels; missing responses stay NaN.
s=lower(strtrim(string(values))); s=erase(s,'num_'); s=erase(s,'key_'); s=erase(s,'button_');
s=regexprep(s,'^[''" ]+|[''" ]+$',''); x=str2double(s);
for i=find(~isfinite(x))'
    token=regexp(char(s(i)),'[-+]?\d+(\.\d+)?','match','once');
    if ~isempty(token), x(i)=str2double(token); end
end
x(ismember(s,["yes","y","true","positive"]))=2;
x(ismember(s,["no","n","false","negative"]))=1;
end
