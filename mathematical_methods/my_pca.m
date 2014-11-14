function [x_bar, P_x, B] = my_pca(X, thresh)
[N d] = size(X);
x_bar = mean(X);
D = (X-repmat(x_bar, N, 1))';
T = D' * D / (N - 1);
S = D * D' / (N - 1);
[U V] = eig(T);
%[U V] = eig(S);

%sort the eigen vectors and values
sort_eig = sortrows([diag(V), (D*U)'], -1);

eig_vals = sort_eig(:,1);
eig_vecs = (sort_eig(:, 2:end))';

%compute t, such that the first t eigen values account thresh% of the
%variance
t = find(cumsum(eig_vals/sum(eig_vals)) >= thresh, 1);
eig_vals = eig_vals(1:t);
P_x = eig_vecs(:, 1:t)./ repmat(sum(eig_vecs(:, 1:t).^2), d, 1);

B = P_x' * D;

