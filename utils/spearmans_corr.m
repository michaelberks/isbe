function rho = spearmans_corr(var1, var2)

if (nargin==0 && nargout==0), test_script(); return; end

ranks1 = rank_vars(var1);
ranks2 = rank_vars(var2);

rho = corr_coeff(ranks1, ranks2);


function ranks = rank_vars(values)

% Get rank of every value in the list
[ignore, tmp] = sort(values);
[ignore, ranks] = sort(tmp);

% Find duplicate ranks and replace with mean rank
unique_values = unique(values);
for v = 1:length(unique_values)
    ranks(values==unique_values(v)) = mean(ranks(values==unique_values(v)));
end


function cc = corr_coeff(var1, var2)

nvar1 = var1 - mean(var1);
nvar2 = var2 - mean(var2);

cc = nvar1(:)'*nvar2(:) / ...
     sqrt(nvar1(:)'*nvar1(:) * nvar2(:)'*nvar2(:));


%% Test script
function test_script()
clc;

v1 = 1:10;
v2 = 1:10;

disp(spearman_corr(v1, v2));

v2 = 10:-1:1;

disp(spearman_corr(v1, v2));

v1 = linspace(-3,3,101);
v2 = 1 ./ (1 + exp(-v1));

disp(spearman_corr(v1, v2));

