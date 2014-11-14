function [cc_mean_shape cc_shapes_a bad_shapes] = compute_cc_mean(cc_shapes, cc_areas, debug_mode)
%COMPUTE_CC_MEAN *Insert a one line summary here*
%   [cc_mean_shape] = compute_cc_mean(cc_shapes)
%
% Inputs:
%      cc_shapes- *Insert description of input variable here*
%
%
% Outputs:
%      cc_mean_shape- *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 17-Nov-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester
%set default arguments
if nargin < 3
    debug_mode = 0;
end

%Procrustes align the CC shapes
[cc_shapes_a, cc_scales, cc_mean, cc_rots, cc_trans] = ...
    align_shapes(cc_shapes, 'area', mean(cc_areas), 'shiftOrigin', 0);

[n_shapes n_pts] = size(cc_shapes);
n_pts = n_pts / 2;

% Need to correct rotation so that chest wall is vertical of left-hand side
chest_vec = cc_mean(end,:) - cc_mean(1+2*n_pts/3,:);

%Compute angle of chest wall vector
theta = acos(chest_vec(1) / sqrt(sum(chest_vec.^2)));

%Build rotation matrix to rotate shapes as desired
rho = 3*pi/2 - theta;
rho_m = [cos(rho) -sin(rho); sin(rho) cos(rho)];

%Rotate the mean CC shape
cc_mean_rho = (rho_m * cc_mean')';

% Update the aligned shapes and transformations given the correctly rotated
% mean shape 
for ii = 1:size(cc_shapes_a,1);
    
    cc_shape = reshape(cc_shapes(ii,:), [], 2);
    [dd Z t] = mb_procrustes(cc_mean_rho, cc_shape);
    
    %Z = t.b * Y * t.T + t.c -> Y = (Z + t.c)* t.Ti / t.b 
    cc_shapes_a(ii,:) = [Z(:,1)' Z(:,2)'];
    cc_rots(:,:,ii) = t.T;
    cc_scales(ii) = t.b;
    cc_trans(ii,:) = t.c(1,:);
    
%     if debug_mode
%         cc_shape_a2 = cc_scales(ii)*cc_shape*cc_rots(:,:,ii) + repmat(cc_trans(ii,:), n_pts,1);
%         cc_shape2 = (cc_shape_a2 - repmat(cc_trans(ii,:), n_pts,1)) * inv(cc_rots(:,:,ii)) / cc_scales(ii);
%         
%         figure; 
%         subplot(1,2,1); axis ij equal; hold on;
%         plot(cc_shape(:,1), cc_shape(:,2), 'b');
%         plot(cc_shape2(:,1), cc_shape2(:,2), 'r:');
% 
%         subplot(1,2,2); axis ij equal; hold on;
%         plot(Z(:,1), Z(:,2), 'b');
%         plot(cc_shape_a2(:,1), cc_shape_a2(:,2), 'r:');
%     end
end

%Now we need to check if we have any bad segmentations - we can do this
%automatically by building a PCA shape model and seeing which breast shapes
%are not well predicted by the model

% Build shape model and look at reconstructions
[cc_mean_pca, P_cc, B_cc] = pca(cc_shapes_a, 1);
recon_erros = zeros(n_shapes,1);
for ii = 1:n_shapes
    
    %Reconstruction from 1st principal mode
    cc_shape_new = cc_mean_pca' + (P_cc(:,1)*B_cc(1,ii));
    cc_shape_new = reshape(cc_shape_new, [], 2);
    
    %Original shape
    cc_shape_a = reshape(cc_shapes_a(ii,:), [], 2);
    recon_erros(ii) = sqrt(mean(sum((cc_shape_new-cc_shape_a).^2, 2)));
end

%Define bad shapes as those that lie more than 3 standard deviations from
%the mean reconstruction error
bad_shapes = recon_erros > (mean(recon_erros) + 3*std(recon_erros));
%
%Throw away the bad shapes
cc_shapes_a(bad_shapes,:) = []; 

%Compute the final mean
cc_mean_shape = reshape(mean(cc_shapes_a), [], 2);

if debug_mode
    figure; hold on; axis ij equal;
    for ii = 1:size(cc_shapes_a,1);
        cc_shape_a = reshape(cc_shapes_a(ii,:), [], 2);
        plot(cc_shape_a(:,1), cc_shape_a(:,2));
    end
    plot(cc_mean_shape(:,1), cc_mean_shape(:,2), 'r');
end
    
%
%Some figure we might want if we're debugging
% figure; plot(cc_mean(1:50,1), cc_mean(1:50,2), 'rx'); axis ij equal; hold on;
% plot(cc_mean(51:100,1), cc_mean(51:100,2), 'gx');
% plot(cc_mean(101:150,1), cc_mean(101:150,2), 'yx');
% 
% figure; plot(cc_mean1(1:50,1), cc_mean1(1:50,2), 'rx'); axis ij equal; hold on;
% plot(cc_mean_rho(51:100,1), cc_mean_rho(51:100,2), 'gx');
% plot(cc_mean_rho(101:150,1), cc_mean_rho(101:150,2), 'yx');
% 
% % Plot all aligned shapes on a single axis
    
