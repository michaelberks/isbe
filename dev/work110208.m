%set up sampling arguments for forest
forest = u_load('C:\isbe\asymmetry_project\data\line_orientation_rfs\243311\random_forest.mat');
sampling_args = u_load('C:\isbe\asymmetry_project\data\line_orientation_rfs\243311\sampling_args.mat');
sampling_args_c.num_levels = sampling_args.num_levels;
sampling_args_c.feature_shape = sampling_args.feature_shape;
sampling_args_c.feature_type = sampling_args.feature_type;
sampling_args_c.do_max = sampling_args.do_max;
sampling_args_c.rotate = sampling_args.rotate;
sampling_args_c.win_size = sampling_args.win_size;
sampling_args_c.use_nag = sampling_args.use_nag;
%

%pre-allocate error spaces
n_angles = 36;
n_cons = 8;
mean_ori_errors_rf = zeros(n_angles, n_cons);
mean_ori_errors_g2 = zeros(n_angles, n_cons);

n_backgrounds = 100;
bg_dir = 'C:\isbe\asymmetry_project\data\synthetic_backgrounds\real512\test\bg';

%create circle mask
xy = repmat(-128:127, 256, 1);
circle_mask = xy.^2 + xy'.^2 < 100^2;
%
%loop through orientations
for ori = 1:n_angles
    theta = 180*ori/n_angles;
    
    %loop through cons
    for con = 1:n_cons
        
        display(['Generating errors for ori = ' num2str(theta) '; con = ' num2str(con);]);
        
        %generate line
        [bar_image, label, label_centre, label_orientation] = create_sin_bar(...
            3, 2*con, theta, 256, 256, 0.5, 128, 128);
        
        %Create image mask
        mask = circle_mask & label_centre;
        n_pts = sum(mask(:));
        
        %pre-allocate errors
        line_errors_rf = zeros(n_backgrounds*n_pts,1);
        line_errors_g2 = zeros(n_backgrounds*n_pts,1);
        
        %loop through backgorunds
        for bb = 1:2:n_backgrounds
            
            %load background
            bg_real = u_load([bg_dir zerostr(bb,5) '.mat']);
            
            %Add bar to background
            test_image = bar_image + bg_real(1:256, 1:256);
            
            %********************** RF ***********************************
            
            %Classify test image
            [orientation_rf] = classify_image(...
                'image_in', test_image, ...
                'forest', forest,...
                'mask', mask,...
                'max_size', 256,...
                'sampling_args', sampling_args_c,...
                'forest_type', 'orientation',...
                'decomp_type', 'dt');
            
            %Convert orientations to degree form
            orientation_rf = mod(180*angle(orientation_rf)/pi,180);
            
            %Compute errors and conctenate
            line_errors_rf((bb-1)*n_pts+1:bb*n_pts) = ...
                abs(mb_mod(orientation_rf(mask) - theta, 180));
            
            %********************** G2 ***********************************
            
            %Apply karssemeijer orientation estimation method
            [line_orientation] = karssemeijer_line_detection(test_image);
            
            %Compute errors and conctenate
            line_errors_g2((bb-1)*n_pts+1:bb*n_pts) = ...
                abs(mb_mod(line_orientation(mask) - theta, 180));
                
        end
        
        %Compute mean errors and save
        mean_ori_errors_rf(ori,con) = mean(line_errors_rf);
        save('C:\isbe\asymmetry_project\data\misc\ori_errors_243311.mat', 'mean_ori_errors_rf');
        
        mean_ori_errors_g2(ori,con) = mean(line_errors_g2);
        save('C:\isbe\asymmetry_project\data\misc\ori_errors_g2.mat', 'mean_ori_errors_g2');
        
    end
end
%%
forest_m = u_load('C:\isbe\asymmetry_project\data\line_orientation_rfs\243629\random_forest.mat');
forest_p = u_load('C:\isbe\asymmetry_project\data\line_orientation_rfs\243630\random_forest.mat');

sm = u_load('C:\isbe\asymmetry_project\data\line_orientation_rfs\243629\sampling_args.mat');
sp = u_load('C:\isbe\asymmetry_project\data\line_orientation_rfs\243630\sampling_args.mat');

sampling_args_m.num_levels = sm.num_levels;
sampling_args_m.feature_shape = sm.feature_shape;
sampling_args_m.feature_type = sm.feature_type;
sampling_args_m.do_max = sm.do_max;
sampling_args_m.rotate = sm.rotate;
sampling_args_m.win_size = sm.win_size;
sampling_args_m.use_nag = sm.use_nag;

sampling_args_p.num_levels = sp.num_levels;
sampling_args_p.feature_shape = sp.feature_shape;
sampling_args_p.feature_type = sp.feature_type;
sampling_args_p.do_max = sp.do_max;
sampling_args_p.rotate = sp.rotate;
sampling_args_p.win_size = sp.win_size;
sampling_args_p.use_nag = sp.use_nag;
%

%pre-allocate error spaces
n_angles = 36;
n_cons = 8;
mean_ori_errors_rf_m = zeros(n_angles, n_cons);
mean_ori_errors_rf_p = zeros(n_angles, n_cons);

n_backgrounds = 100;
bg_dir = 'C:\isbe\asymmetry_project\data\synthetic_backgrounds\real512\test\bg';

%create circle mask
xy = repmat(-128:127, 256, 1);
circle_mask = xy.^2 + xy'.^2 < 100^2;
%%
%loop through orientations
for ori = 1:n_angles
    theta = 180*ori/n_angles;
    
    %loop through cons
    for con = [1:3 5:n_cons]
        
        display(['Generating errors for ori = ' num2str(theta) '; con = ' num2str(con);]);
        
        %generate line
        [bar_image, label, label_centre, label_orientation] = create_sin_bar(...
            3, 2*con, theta, 256, 256, 0.5, 128, 128);
        
        %Create image mask
        mask = circle_mask & label_centre;
        n_pts = sum(mask(:));
        
        %pre-allocate errors
        line_errors_rf_p = zeros(n_backgrounds*n_pts,1);
        line_errors_rf_m = zeros(n_backgrounds*n_pts,1);
        
        %loop through backgorunds
        for bb = 1:2:n_backgrounds
            
            %load background
            bg_real = u_load([bg_dir zerostr(bb,5) '.mat']);
            
            %Add bar to background
            test_image = bar_image + bg_real(1:256, 1:256);
            
            %********************** RF MAG ********************************
            
            %Classify test image
            [orientation_rf_m] = classify_image(...
                'image_in', test_image, ...
                'forest', forest_m,...
                'mask', mask,...
                'max_size', 256,...
                'sampling_args', sampling_args_m,...
                'forest_type', 'orientation',...
                'decomp_type', 'dt');
            
            %Convert orientations to degree form
            orientation_rf_m = mod(180*angle(orientation_rf_m)/pi,180);
            
            %Compute errors and conctenate
            line_errors_rf_m((bb-1)*n_pts+1:bb*n_pts) = ...
                abs(mb_mod(orientation_rf_m(mask) - theta, 180));
            
            %********************** RF PHASE ******************************
            
            %Classify test image
            [orientation_rf_p] = classify_image(...
                'image_in', test_image, ...
                'forest', forest_p,...
                'mask', mask,...
                'max_size', 256,...
                'sampling_args', sampling_args_p,...
                'forest_type', 'orientation',...
                'decomp_type', 'dt');
            
            %Convert orientations to degree form
            orientation_rf_p = mod(180*angle(orientation_rf_p)/pi,180);
            
            %Compute errors and conctenate
            line_errors_rf_p((bb-1)*n_pts+1:bb*n_pts) = ...
                abs(mb_mod(orientation_rf_p(mask) - theta, 180));
                
        end
        
        %Compute mean errors and save
        mean_ori_errors_rf_m(ori,con) = mean(line_errors_rf_m);
        mean_ori_errors_rf_p(ori,con) = mean(line_errors_rf_p);
        
        save('C:\isbe\asymmetry_project\data\misc\ori_errors_243629.mat', 'mean_ori_errors_rf_m');
        save('C:\isbe\asymmetry_project\data\misc\ori_errors_243630.mat', 'mean_ori_errors_rf_p');
        
    end
end