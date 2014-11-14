function [mean_X, P_x, B_x, L_x,...
    mean_scale, P_scale, B_scale, L_scale mean_target, shape_scale, X_pro]...
    = shape_model(X, area, thresh, seed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% author:   Michael Berks
% date:     26/05/2006  10:30
%
% function: take input of concatenate shape vectors and computes shape model using PCA
%           each row of X is shape vector of form: (x1,x2,...,xn,y1,y2,...,yn)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[N dim] = size(X); dim = dim/2; % N is number of shape vectors, dim is number of points in each vector
if nargin < 4;
    seed = 1; %index of intial shape used to align
end
%
% first align vectors using procrustes analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X_pro = zeros(N, 2*dim);
shape_scale = zeros(N, 1);
sq_diff = zeros(N,1);

for ii=1:N
    shape_vec(:,1) = X(ii, 1:dim);
    shape_vec(:,2) = X(ii, dim+1:end);
    [dd Z t] = procrustes([X(seed,1:dim)' X(seed,dim+1:end)'], shape_vec);
    X_pro(ii, :) = [Z(:,1)', Z(:,2)'];
    shape_scale(ii) = t.b;
    sq_diff(ii) = dd;
end

mean_sq_diff = mean(sq_diff);
display(['mean change = ', num2str(mean_sq_diff)]);

go_on = 1; jj = 1;
while go_on
    
    mean_shape = mean(X_pro);
    mean_shape = reshape(mean_shape, dim, 2);
    mean_shape = mean_shape - repmat(mean(mean_shape), dim, 1);
    scale_factor = sqrt(area/polyarea(mean_shape(:,1), mean_shape(:,2)));
    mean_shape = mean_shape*scale_factor;

    th = atan2(mean_shape(1,1), mean_shape(1,2));
    mean_shape = ([cos(th) -sin(th); sin(th) cos(th)] * mean_shape')';
    
    for ii = 1:N
        [dd Z t] = procrustes(mean_shape,...
            [X_pro(ii, 1:dim)' X_pro(ii, dim+1:end)']);
        X_pro(ii, :) = [Z(:,1)', Z(:,2)'];
        sq_diff(ii) = dd;
        shape_scale(ii) = shape_scale(ii)*t.b;
        clear dd Z t;
    end
    go_on = abs(mean_sq_diff - mean(sq_diff)) > 1e-4;
    mean_sq_diff = mean(sq_diff);
    display(['mean change = ', num2str(mean_sq_diff)]);
    
    jj = jj+1;
end
mean_target = mean_shape;
display(['stopped after ', num2str(jj-1), ' iterations']);

%
% the perform PCA on the aligned points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[mean_X, P_x, B_x, L_x] = pca(X_pro, thresh);%, 0);
[mean_scale, P_scale, B_scale, L_scale] = pca(shape_scale, thresh);%, 0);

