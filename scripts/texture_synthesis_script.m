% function texture_synthesis_script(id_idx, data_type)
% %final cluster wrapper script
% 
% %id_idx = 1,...,30 and comes from $SGE_IDX batch job variable
% [levels window_sizes] = meshgrid(2:5, [5 7 9 11]);
% 
% 
% pyramid = u_load([mberksroot, 'background/pyramid/normal512/o04_001LCC_1024_3427_865_pyramid']);
% 
% syn_args.PathToTextureGMM = [mberksroot, 'background/models/window_size/',...
%     data_type, '_level_',...
%     num2str(levels(id_idx)) '_win_',...
%     num2str(window_sizes(id_idx))];
% 
% syn_args.TargetImage = pyramid{levels(id_idx), 1};
% syn_args.FilledImage = logical(ones(2^(11-levels(id_idx)))); %#ok
% 
% from = 1 + 2^(9 - levels(id_idx));
% to = (2^(11-levels(id_idx))) - (2^(9 - levels(id_idx)));
% 
% syn_args.FilledImage(from:to, from:to) = 0;
% 
% syn_args.SaveFile = [mberksroot, 'background/syn/window_size/',...
%     'synthesis_level_',...
%     num2str(levels(id_idx)) '_win_',...
%     num2str(window_sizes(id_idx))];
% 
% syn_args.MovieFilename = 'C:\isbe\dev\misc\test_gif.gif';
% 
% display(syn_args);
% syn_image = mb_gmm_tex_synthesis(syn_args);
% 
% write_im_from_colormap(syn_image, [syn_args.SaveFile(1:end-4), '.bmp'], jet(256), [-4 4]);
% clear;
% exit
% 
% pyr_args.TargetPyramid = u_load('C:\isbe\dev\background\pyramid\normal_2\normal001.bmp_pyramid.mat');
% pyr_args.FilledImage = logical(ones(size(pyr_args.TargetPyramid{2,1})));
% row_centre = round(size(pyr_args.TargetPyramid{1}, 1) / 2);
% col_centre = round(size(pyr_args.TargetPyramid{1}, 2) / 2);
% pyr_args.FilledImage(row_centre-64:row_centre+64-1, col_centre-64:col_centre+64-1) = 0;
% 
% pyr_args.ModelDir = 'C:\isbe\dev\background\results\2levels\';
% pyr_args.ModelName = 'normal_2_model';
% 
% [synthesised_image, pyramid, cluster_image] = mb_gmm_pyr_synthesis(pyr_args);
%%
display(['--Texture synthesis script started: ' datestr(now)]);
clear

for ii = 7
    pyr_args.TargetPyramid = u_load([mberksroot,...
        'background/old_pyramid/normal_2/normal', zerostr(ii, 3), '.bmp_pyramid.mat']);
    [rows cols] = size(pyr_args.TargetPyramid{2,1});
    row_centre = round(rows / 2);
    col_centre = round(cols / 2);
%     pyr_args.FilledImage = logical(ones(rows, cols)); %#ok
%     pyr_args.FilledImage(row_centre-64:row_centre+63, col_centre-64:col_centre+63) = 0;

    % Make a the biggest circular mask
    m = min([rows cols]);
    rad = floor((m - 128) / 2);
    [x y] = meshgrid(1:cols, 1:rows);
    pyr_args.FilledImage = (x - col_centre).^2 + (y - row_centre).^2 > rad.^2;
    clear x y m rad row_centre col_centre;

    pyr_args.ModelName = 'normal_2_k10_c_model';

    pyr_args.ModelDir = [mberksroot, 'background/models5/normal_2_2levels/'];
    pyr_args.CutOffLevel = 4;
    pyr_args.ConditionLevels = 1;
    pyr_args.SaveFile = [mberksroot, 'background/syn/epsrc/conditioned_synthesis_circle_6_', zerostr(ii, 3)];
    pyr_args.MovieFilename = [mberksroot, 'background/syn/epsrc/conditioned_synthesis_circle_6_', zerostr(ii, 3)];
    mb_gmm_pyr_synthesis(pyr_args);

    clear pyr_args;
end

% clear
display(['--Texture synthesis script finished: ' datestr(now)]);
% exit