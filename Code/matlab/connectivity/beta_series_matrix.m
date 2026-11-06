function [R,Z,bad_roi] = beta_series_matrix(X)
% BETA_SERIES_MATRIX Pearson correlation and Fisher-z of trial-wise ROI means.
% Constant/nonfinite ROI series produce unavailable edges. The correlation
% is clipped to +/-0.999999 before atanh; stored matrix diagonals are zero.
validateattributes(X,{'numeric'},{'2d','real','nonempty'});
X = double(X);
bad_roi = std(X,0,1)<eps | any(~isfinite(X),1);
X(:,bad_roi) = NaN;
R = corr(X,'Rows','pairwise','Type','Pearson');
R(1:size(R,1)+1:end) = 0;
Z = atanh(min(max(R,-0.999999),0.999999));
Z(1:size(Z,1)+1:end) = 0;
end
