function [x_bar, P_x, B, eig_vals] = pca(X, thresh, plot_modes)

if nargin < 3
    plot_modes = 0;
end

[N d] = size(X);
x_bar = mean(X);
D = (X-repmat(x_bar, N, 1))';

%S = D * D' / (N - 1);
%[U V] = eig(S);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = D' * D / (N); %using T is more effcient for N < d
[U V] = eig(T);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%sort the eigen vectors and values
sort_eig = sortrows([diag(V), U'], -1);

eig_vals = sort_eig(:,1);
eig_vecs = (sort_eig(:, 2:end))';

%compute t, such that the first t eigen values account thresh% of the
%variance

if thresh > 1
    t = round(thresh);
elseif thresh > 0
    t = find(cumsum(eig_vals) >= sum(eig_vals)*thresh, 1);

    if isempty(t)
       t = 1;
    elseif(t == N)
       t = N - 1;
    end

else
    t = N - 1;
    Warning('Threshold should be between 0 and 1: all modes kept');
end

eig_vals = eig_vals(1:t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_x = (D*eig_vecs(:, 1:t))./ repmat(sqrt(N*eig_vals'), d, 1);
% use if use T as cov mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%P_x = eig_vecs(:, 1:t)./ repmat(sum(eig_vecs(:, 1:t).^2), d, 1);
B = P_x' * D;

if plot_modes
    
    mode1_plus_2sd = x_bar + (2*P_x(:,1)*sqrt(eig_vals(1)))';
    mode1_min_2sd = x_bar + (-2*P_x(:,1)*sqrt(eig_vals(1)))';
    
    mode2_plus_2sd = x_bar + (2*P_x(:,2)*sqrt(eig_vals(2)))';
    mode2_min_2sd = x_bar + (-2*P_x(:,2)*sqrt(eig_vals(2)))';
    
    figure,
    subplot(2,3,1)
    plot(mode1_min_2sd(1:end/2), mode1_min_2sd(1+end/2:end)); 
    subplot(2,3,2)
    plot(x_bar(1:end/2), x_bar(1+end/2:end)); 
    subplot(2,3,3)
    plot(mode1_plus_2sd(1:end/2), mode1_plus_2sd(1+end/2:end));
    
    subplot(2,3,4)
    plot(mode2_min_2sd(1:end/2), mode2_min_2sd(1+end/2:end));
    subplot(2,3,5)
    plot(x_bar(1:end/2), x_bar(1+end/2:end)); 
    subplot(2,3,6)
    plot(mode2_plus_2sd(1:end/2), mode2_plus_2sd(1+end/2:end));
    figure,
    bar(eig_vals);
end
    
    