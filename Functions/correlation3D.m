function r = correlation3D(a,b)

% Synthetic data for testing
b(1, 1, :) = 3 * a(1, 1, :) - 2; % r(1,1) should be + 1;
b(1, 2, :) = -17 * a(1, 2, :) + 8; % r(1,2) should be - 1;
% rest of r should be random between +1 and -1
% Compute correlations on third dimension
% Remove means 
az = bsxfun(@minus, a, mean(a,3));
bz = bsxfun(@minus, b, mean(b,3));
% Standard Pearson correlation coefficient formula
a2 = az .^ 2;
b2 = bz .^ 2;
ab = az .* bz;
r = sum(ab, 3) ./ sqrt(sum(a2, 3) .* sum(b2, 3));