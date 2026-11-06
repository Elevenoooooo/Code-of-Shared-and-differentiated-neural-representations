function [D, trialKeep, nVoxel] = rsa_neural_rdm(patterns, minTrials, minVoxels)
%RSA_NEURAL_RDM One minus Pearson correlation of trial-wise voxel patterns.
% Rows are LSS trials and columns are voxels. Finite, nonconstant voxels
% contribute to the correlation calculation within each ROI.
if nargin < 2, minTrials = 20; end
if nargin < 3, minVoxels = 10; end
patterns = double(patterns); trialKeep = any(isfinite(patterns), 2);
assert(sum(trialKeep) >= minTrials, 'Fewer than %d valid trials.', minTrials);
X = patterns(trialKeep, :);
sd = std(X, 0, 1, 'omitnan');
keep = all(isfinite(X), 1) & isfinite(sd) & sd > eps;
nVoxel = sum(keep);
assert(nVoxel >= minVoxels, 'Fewer than %d finite nonconstant voxels.', minVoxels);
X = X(:, keep);
D = 1 - corr(X', 'Rows', 'pairwise', 'Type', 'Pearson');
D = (D + D') / 2; D(1:size(D,1)+1:end) = 0;
end
