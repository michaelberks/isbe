function [] = pyr_synthesis(start_idx, data_type)

pyr_list = dir([mberksroot, 'background/pyramid/' data_type , '/*pyr*']);


display(['--Texture synthesis script started: ' datestr(now)]);

for ii = start_idx:start_idx+19
    
    if ii > length(pyr_list); break; end;
    
    pyr_args.TargetPyramid = u_load([mberksroot,...
        'background/pyramid/' data_type , '/', pyr_list(ii).name]);
    
    [rows cols] = size(pyr_args.TargetPyramid{2,1});
    row_centre = round(rows / 2);
    col_centre = round(cols / 2);

    % Make a the biggest circular mask
    m = min([rows cols]);
    rad = floor((m - 128) / 2);
    [x y] = meshgrid(1:cols, 1:rows);
    pyr_args.FilledImage = (x - col_centre).^2 + (y - row_centre).^2 > rad.^2;
    clear x y m rad row_centre col_centre;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pyr_args.ModelName = 'normal512_k10_w1_11_w2_5';
    pyr_args.ModelDir = [mberksroot, 'background/models/pyramid/'];
    pyr_args.CutOffLevel = 5;
    pyr_args.ConditionLevels = 1;
    pyr_args.SaveFile = [mberksroot, 'background/syn/pyramid/',...
        data_type, '_', pyr_args.ModelName, '_' zerostr(ii, 3)];
    mb_gmm_pyr_synthesis(pyr_args);

    clear pyr_args;
end

% clear
display(['--Texture synthesis script finished: ' datestr(now)]);