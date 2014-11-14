function [] = test_orientation_images4(ii, ff, varargin)
%TEST_ORIENTATION_IMAGES *Insert a one line summary here*
%   [] = test_orientation_images(ii, ff)
%
% Inputs:
%      ii - *Insert description of input variable here*
%
%      ff - *Insert description of input variable here*
%
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 11-Jan-2011
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
figure(ff);

warning('off', 'load_uint8:missing_variables');
test_dir = [asymmetryroot 'data\synthetic_lines\real512\'];
prob_dir = [asymmetryroot 'data\synthetic_lines\real512\results\'];
param_dirs = {'233902', '286712', 'g2d_scales_16\orientations\'};
titles = {'RF/DT-CWT, real backgrounds', 'RF/DT-CWT, synthetic backgrounds', 'Gaussian 2nd derivatives'};

spacing = 1;
num_angles = 36;
ang_res = 180 / num_angles;
arrow_colors = hsv(num_angles);


test_im = u_load([test_dir 'image' zerostr(ii,3) '.mat']);
label = load([test_dir 'labels\label' zerostr(ii,3) '.mat']);

a = zeros(4,1);

for kk = 1:4
    if kk == 4
        orientation_map = label.label_orientation;
        a(kk) = subplot(2,2,kk); imagesc(label.label); axis image; colormap(gray(256)); hold on;
        mask = label.label > 0;
        
    else
        orientation_map = load_uint8([prob_dir param_dirs{kk} '\image' zerostr(ii,3) '_class.mat']);
        if ~isreal(orientation_map)
            orientation_map = mod(180*angle(orientation_map)/pi,180);
        end
        a(kk) = subplot(2,2,kk); imagesc(test_im); axis image; colormap(gray(256)); hold on;
        mask = label.label_centre > 0;
        
        error_mask = mask & label.label < 2; 
        ori_errors = mb_mod(orientation_map(error_mask) - label.label_orientation(error_mask), 180);
        title(titles{kk});
        xlabel(['Mean orientation error = ' num2str(mean(abs(ori_errors)))]);
    end
    
    
    [r c] = size(mask);

    spacing_mask = false(r, c);
    spacing_mask(1:spacing:r, 1:spacing:c) = true;
    mask = mask & spacing_mask;

    for jj = 0:num_angles
        theta = (jj - 0.5)*ang_res;

        %Get mask of pixels that have orientation within theta range
        angle_mask = mask & ...
             (orientation_map > theta - 0.5*ang_res) &...
             (orientation_map <= theta + 0.5*ang_res);

        [y x] = find(angle_mask);
        u = cos(pi*orientation_map(angle_mask)/180);
        v = -sin(pi*orientation_map(angle_mask)/180);

        quiver(a(kk), x, y, 4*u, 4*v, 0, 'color', arrow_colors(mod(jj, num_angles)+1,:));
    end
end
linkaxes(a);
zoom on;

hManager = uigetmodemanager(ff);
set(hManager.WindowListenerHandles,'Enable','off');
