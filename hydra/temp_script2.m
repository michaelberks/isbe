function temp_script2(idx, num_jobs)


%--------------------------------------------------------------------------

display(['temp_script2 running: ' datestr(now)]);

if ~ispc
     [z num_jobs] = unix('echo $NUM_JOBS'); num_jobs = str2num(num_jobs);
end

ret_dir = [asymmetryroot 'data/retinograms/DRIVE/'];

num_images = 40;
images_per_job = ceil(num_images / num_jobs);
start_idx = (idx-1)*images_per_job + 1;
end_idx = min(num_images, idx*images_per_job);

sampling_args.feature_shape = 'rect';
sampling_args.feature_type = 'complex';
sampling_args.levels = 1:6;

for ii = start_idx:end_idx
    display(['processing image ' num2str(ii)]);
    if ii < 21
        data = 'test';
    else
        data = 'training';
    end
    
    %Load in ground truth
    gt = logical(imread([ret_dir data '\1st_manual\' zerostr(ii,2) '_manual1.gif']));
    
    %Load in retinogram
    ret = imread([ret_dir data '\images\' zerostr(ii,2) '_' data '.tif']);
    
    %load in mask of edge of image
    mask = logical(imread([ret_dir data '\mask\' zerostr(ii,2) '_' data '_mask.gif']));
    mask_inner = imerode(mask, strel('disk', 10));
    mask_inner2 = imerode(mask, strel('disk', 20));
    mask_ring = mask_inner & ~mask_inner2;
    
    %Load in the orientation mask
    ori_map = u_load([ret_dir data '\orientations\' zerostr(ii,2) '_ori1.mat']);
    
    %Get row and columns of pixels in ground truth
    [rows cols] = find(gt);
    
    %Sample orientations and save
    y = ori_map(gt);
    save([asymmetryroot 'data/synthetic_data/drive/dt/' data '/y_' zerostr(ii,2) '.mat'], 'y');
    
    clear y;
    
    %Get positions to extend border
    [y_int x_int] = find(~mask_inner);
    [y_o x_o] = find(mask_ring);
    
    %Work on each RGB channel idenpendently
    rgb = 'rgb';
    for ch = 1:3
        channel = ret(:,:,ch);
        
        %Process the retinogram to get rid of the hard edge
        z_o = double(channel(mask_ring));
        z_int = griddata(x_o, y_o, z_o, x_int, y_int, 'nearest');
        channel(~mask_inner) = z_int;
        clear z_int;
        
        %Compute dual-tree of retinogram
        dt = compute_dual_tree(channel, 6, 0);
    
        %Sample DT coefficients from specified rows and cols according to
        %sampling arguments
        X = sample_dt_data(dt, rows, cols, sampling_args);
        X = reshape(X, size(X,1), 9, 6, 6);

        %Save X
        save([asymmetryroot 'data/synthetic_data/drive/dt/' data '/X_' rgb(ch) '_' zerostr(ii,2) '.mat'], 'X');
        clear X;
    end  
end
display(['temp_script2 completed: ' datestr(now)]);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Old stuff
%--------------------------------------------------------------------------
%
% im = double(imread([asymmetryroot 'data/misc/o04_010RCC_1024_3797_3365.bmp']));
% xx = repmat(1:512, 512, 1);
% mask = (xx-256).^2 + (xx' - 256).^2 < 256^2;
% %for theta = 9*(0:19)
% 
%     theta = idx;%9*(idx-1);
% 
%     i1 = imrotate(im, theta, 'crop');
%     i1 = i1(257:768, 257:768);
% 
%     i2 = imrotate(im, theta+1, 'crop');
%     i2 = i2(257:768, 257:768);
% 
%     args.image1 = i1;
%     args.image2 = i2;
%     args.mask1 = mask;
%     args.mask2 = mask;
%     args.num_train = 10000;
%     args.num_levels = 5;
%     args.tree_dir = [asymmetryroot 'data/misc/'];
%     args.do_test1 = 1;
%     args.do_test2 = 1;
%     args.use_probs = 0;
%     args.feature_type = 'ilp';
%     args.do_max = 0;
%     args.win_size = 1;
%     args.max_test_size = 64;
%     args.n_trees = 100;
%     args.d = round(args.win_size * sqrt(2*args.num_levels * (6 - 5*args.do_max)));
%     
%     args.save_path = [asymmetryroot 'data/misc/same_image_rf_W' num2str(args.win_size) 'L' num2str(args.num_levels) 'M' num2str(args.do_max) 't' zerostr(2*theta, 3) '.mat'];
%     clear i1 i2 im xx mask;
%     %
%     mb_random_forest_two_images(args);
% 
% %end