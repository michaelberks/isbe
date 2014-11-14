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
% %%
display(['--Texture synthesis script started: ' datestr(now)]);
clear

for ii = 1:15
    pyr_args.TargetPyramid = u_load([mberksroot,...
        'background/g_pyramid/mass_2/mass', zerostr(ii, 3), '_pyramid.mat']);
    [rows cols] = size(pyr_args.TargetPyramid{1,1});
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

    pyr_args.ModelName = 'mass_2_g_20_model';
    pyr_args.ModelNameCls = 'mass_2_g_cls_20_model';
    cls_result = u_load(['C:/isbe/dev/background/cls/mass_2/mass_cls_', zerostr(ii, 3), '.mat']);
    pyr_args.ClsMap = mb_cls_map_merge(cls_result);
    
    pyr_args.ModelDir = [mberksroot, 'background/results/g_pyramid/'];
    pyr_args.CutOffLevel = 4;
    pyr_args.CutOffLevelCls = 3;
    pyr_args.ConditionLevels = 1;
    pyr_args.SaveFile = [mberksroot, 'background/syn/g_pyramid/gp_circle', zerostr(ii, 3)];

    mb_gmm_g_pyr_synthesis_cls(pyr_args);

    clear pyr_args;
end

clear
display(['--Texture synthesis script finished: ' datestr(now)]);
exit