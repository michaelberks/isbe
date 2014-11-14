data_type = 'mass_2';
dt_list = dir([mberksroot, 'background/dual_tree/' data_type , '/*dual_tree*']);


display(['--Texture synthesis script started: ' datestr(now)]);

for ii = 1:1
    
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
    dt_args.FilledImage = (x - col_centre).^2 + (y - row_centre).^2 > 40.^2;
    clear x y m rad row_centre col_centre;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dt_args.ModelDir = [mberksroot, 'background/models/dual_tree/'];
    dt_args.ModelName = 'mass_2_k10_w1_3_w2_1';
    dt_args.SynthesisMode = 'patch-wise';
    dt_args.MovieFilename = 'test_movie.gif';
    dt_args.SaveFile = [mberksroot, 'background/syn/dual_tree/',...
        data_type, '_', dt_args.ModelName, '_' zerostr(ii, 3)];

    [synthesised_image, dual_tree] = mb_gmm_dual_tree_synthesis_final(dt_args);

    %clear dt_args;
end

display(['--Texture synthesis script finished: ' datestr(now)]);