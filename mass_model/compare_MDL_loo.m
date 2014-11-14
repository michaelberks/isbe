% COMPARE_MDL_LOO put mass outlines into shape matrix 
%    [shapes] = compare_MDL_loo(model_shapes, mass_path, npt)
%
%    inputs:
%       model_shapes  - 
%       mass_path   - File path to mass folder
%       n_pts       - number of points to use per shape          
%                       D is the dimensionality of the data points 
%
%    outputs:
%       
%
%    notes:
%
% Copyright: (C) 2006-2008 University of Manchester
% Author: Michael Berks

function [loo_errors_pro, loo_errors_mdl, shapes_m] =...
    compare_MDL_loo(shapes, size_shape_vec, varargin)

default.mass_path = 'C:\isbe\dev\masses\';
default.mass_files = 'C:\isbe\dev\u_files.mat';
default.n_its = 10;
args = u_packargs(varargin, 0, default);

%Load model shapes and standardise n_pts


if isempty(shapes);
    temp = load(args.mass_files);
    mass_files = temp.u_files1;    
    [shapes_m.unaligned, mass_areas] = get_shapes_from_masses(mass_files,...
        size_shape_vec, 'mass_path', args.mass_path);
else
    shapes_m.unaligned = shapes; clear shapes;
end

N = size(shapes_m.unaligned, 1);    

% Align shapes using procrustes method
%%%
shapes_m.pro = align_shapes(shapes_m.unaligned, 20000);

% Reshape the matrix into 2 x size_shape_vec x N matrix for MDL
shapes_m.unaligned_mdl = zeros(2,size_shape_vec, N);
for ii = 1:N; 
    shapes_m.unaligned_mdl(:,:,ii) = [shapes_m.pro(ii,1:end/2);...
                                      shapes_m.pro(ii,end/2+1:end)];
end

% Align shapes using the MDL model
%%%
[r_pts] = amb_automodelbuild (shapes_m.unaligned_mdl,...
    'circle', 'nIterations', args.n_its, 'saveFrequency', 10,...
    'Quiet', 0, 'optimisePose', 0, 'optimiseOrigin', 1,...
    'nExamplesToOptimise', N,...
    'initialAlign', 0, 'initialOriginOptimisation', 0);
% Reshape the MDL reparametrised points
shapes_m.mdl = zeros(N, 2*size_shape_vec);
for ii = 1:N; shapes_m.mdl(ii,:) = [r_pts(1,:,ii) r_pts(2,:,ii)]; end

m_mdl_shapes = reshape(mean(shapes_m.mdl), [], 2);
m_pro_shapes = reshape(mean(shapes_m.pro), [], 2);

figure; hold on;
plot(m_mdl_shapes(:,1), m_mdl_shapes(:,2));
plot(m_pro_shapes(:,1), m_pro_shapes(:,2), 'r:');
clear m_pro_shapes m_mdl_shapes

% % Align shapes using the MDL model
% %%%
% [r_pts] = amb_automodelbuild (shapes_m.unaligned_mdl,...
%     'circle', 'nIterations', args.n_its, 'saveFrequency', 10,...
%     'Quiet', 0, 'optimisePose', 1, 'optimiseOrigin', 1,...
%     'nExamplesToOptimise', N,...
%     'initialAlign', 1, 'initialOriginOptimisation', 1);
% % Reshape the MDL reparametrised points
% shapes_m.mdl = zeros(N, 2*size_shape_vec);
% for ii = 1:N; shapes_m.mdl(ii,:) = [r_pts(1,:,ii) r_pts(2,:,ii)]; end
% 
% % Align shapes using procrustes method with area set to mean area of MDL
% % shapes
% %%%
% m_mdl_shapes = reshape(mean(shapes_m.mdl), [], 2);
% 
% figure; hold on;
% plot(m_mdl_shapes(:,1), m_mdl_shapes(:,2));
% 
% shapes_m.pro = align_shapes(shapes_m.unaligned,...
%     polyarea(m_mdl_shapes(:,1), m_mdl_shapes(:,2)));
% 
% m_pro_shapes = reshape(mean(shapes_m.pro), [], 2);
% plot(m_pro_shapes(:,1), m_pro_shapes(:,2), 'r:');
% clear m_pro_shapes m_mdl_shapes

%
% Leave-one-out model building
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
loo_errors_pro = zeros(N, 1);
loo_errors_mdl = zeros(N, 1);

for ii = 1:N
    %remove i-th shape from data and set as unseen test shape
    loo_shapes_pro = shapes_m.pro([1:ii-1, ii+1:end], :);
    loo_shapes_mdl = shapes_m.mdl([1:ii-1, ii+1:end], :);
    
    shape_t.pro = shapes_m.pro(ii, :);
    shape_t.mdl = shapes_m.mdl(ii, :);
    
    % Compute model parameters for Pro and MDL shapes
    [mean_pro, P_pro] = pca(loo_shapes_pro, 0.98);
    [mean_mdl, P_mdl] = pca(loo_shapes_mdl, size(P_pro,2));
    
    % Compute shape error for two models
    B_pro = P_pro' * (shape_t.pro - mean_pro)';
    B_mdl = P_mdl' * (shape_t.mdl - mean_mdl)';

    shape_rg.pro = mean_pro + (P_pro*B_pro)';
    shape_rg.mdl = mean_mdl + (P_mdl*B_mdl)';
    
    % Compute loo errors a RMS difference between regenerated shape and
    % test shape
    loo_errors_pro(ii) = sqrt(mean(sum((reshape(shape_rg.pro, [], 2) -...
                            reshape(shape_t.pro, [], 2)).^2, 2)));
    loo_errors_mdl(ii) = sqrt(mean(sum((reshape(shape_rg.mdl, [], 2) -...
                            reshape(shape_t.mdl, [], 2)).^2, 2)));
    
    % MDL values?
end 
    
    
    