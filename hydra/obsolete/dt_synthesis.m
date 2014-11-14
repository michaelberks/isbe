function [] = dt_synthesis(start_idx, data_type, model_name)

dt_list = dir([mberksroot, 'background/dual_tree/' data_type , '/*dual_tree*']);


display(['--Texture synthesis script started: ' datestr(now)]);

for ii = start_idx:start_idx+19
    
    if ii > length(dt_list); break; end;
    
    dt_args.TargetDualTree = u_load([mberksroot,...
        'background/dual_tree/', data_type, '/', dt_list(ii).name]);
    
    [rows cols] = size(dt_args.TargetDualTree{1}(:,:,1));
    row_centre = round(rows / 2);
    col_centre = round(cols / 2);

    % Make a the biggest circular mask
    m = min([rows cols]);
    rad = floor((m - 128) / 2);
    [x y] = meshgrid(1:cols, 1:rows);
    dt_args.FilledImage = (x - col_centre).^2 + (y - row_centre).^2 > rad.^2;
    clear x y m rad row_centre col_centre;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dt_args.ModelDir = [mberksroot, 'background/models/dual_tree/'];
    dt_args.ModelName = model_name;
    dt_args.ConditionLevels = 0;
    dt_args.SaveFile = [mberksroot, 'background/syn/dual_tree/',...
        data_type, '_', dt_args.ModelName, '_' zerostr(ii, 3)];

    mb_gmm_dual_tree_synthesis(dt_args);

    clear dt_args;
end

% clear
display(['--Texture synthesis script finished: ' datestr(now)]);