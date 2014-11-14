function [] = initialise_aam_candidates(image_path, vessel_width, vessel_ori, candidates_xy, varargin)

args = u_packargs(varargin, '0', ...
    'aam_dir', [nailfoldroot 'data/rsa_study/test/aam'],...
    'mean_shape', [],...
    'base_width', 15,...
    'width_sigma', 2,...
    'ori_sigma', 0,...
    'max_num_candidates', inf,...
    'thetas', 0,...linspace(-pi/8, pi/8, 10),...
    'scales', [0.8 0.9 1.0 1.1 1.2],...
    'debug', 0);
clear varargin;

%Load in nailfold
[~,~,file_ext] = fileparts(image_path);
if strcmp(file_ext, '.mat')
    nailfold = u_load(image_path);
else
    nailfold = imread(image_path); 
end
nailfold = nailfold(:,:,1);

if args.width_sigma
    g_width = gaussian_filters_1d(args.width_sigma);
    g_width = g_width / sum(g_width);
    vessel_width = conv2(g_width', g_width, vessel_width, 'same');
end
if args.ori_sigma
    g_ori = gaussian_filters_1d(args.ori_sigma);
    g_ori = g_ori / sum(g_ori);
    vessel_ori = conv2(g_ori', g_ori, vessel_ori, 'same');
end

%--------------------------------------------------------------------------
%Now work out an initial alignment for each maxima and write out this
%candidate
num_candidates = min( size(candidates_xy,1), args.max_num_candidates);

%Make a folder to store the candidates in
create_folder(args.aam_dir);
create_folder([args.aam_dir '\images\']);
create_folder([args.aam_dir '\points\']);
create_folder([args.aam_dir '\out_points\']);
candidates_fid = initialise_candidates_file(args.aam_dir);

num_pts = size(args.mean_shape, 1);

%We need the orientations of the mean shape if we're going to try and
%optimise a local fit
if ~isempty(args.thetas)
    
    [nxy] = compute_spline_normals(args.mean_shape);
    mean_xy_ori = exp(-2i*atan(nxy(:,1)./nxy(:,2)));
    
    %May as well get these too while we're here...
    num_scales = length(args.scales);
end

display(['Matching pose for ' num2str(num_candidates) ' apex candidates']);
for i_ap = 1:num_candidates
    display(['Matching pose for candidate ' num2str(i_ap)]);

    ax = round(candidates_xy(i_ap,1));
    ay = round(candidates_xy(i_ap,2));
    
    theta = angle(vessel_ori(ay, ax)) / 2;
    
    test_thetas = args.thetas;
    if theta > pi/2.5 && ~isempty(args.thetas)
        test_thetas = [args.thetas args.thetas + pi];        
    end
    
    scale = vessel_width(ay, ax) / args.base_width;
    
    rot = [cos(theta) -sin(theta); sin(theta) cos(theta)]; 
    %rot = [1 0; 0 1];
    vessel_xy = args.mean_shape * rot * scale;
    
    %Compute position of RoI about candidate
    sr = max(1, floor(min(vessel_xy(:,2)) - 50 + candidates_xy(i_ap,2)));
    er = min(size(nailfold,1), ceil(max(vessel_xy(:,2)) + 50 + candidates_xy(i_ap,2)));
    sc = max(1, floor(min(vessel_xy(:,1)) - 50 + candidates_xy(i_ap,1)));
    ec = min(size(nailfold,2), floor(max(vessel_xy(:,1)) + 50 + candidates_xy(i_ap,1)));

    %Sample patch from image
    image_patch = nailfold(sr:er, sc:ec);  
    ori_patch = vessel_ori(sr:er, sc:ec);
    
    offset_x = candidates_xy(i_ap,1) - sc;
    offset_y = candidates_xy(i_ap,2) - sr;
    
    if isempty(args.thetas)
        
        %If we're not optimising, just apply the offsets to the scaled and
        %rotated vessel
        vessel_xy(:,1) = vessel_xy(:,1) + offset_x;
        vessel_xy(:,2) = vessel_xy(:,2) + offset_y;
        
    else
    
        %Transform vessels xy orientations
        vessel_xy_ori = mean_xy_ori * exp(-2i*theta);
        
        best_fit = -inf;
        
        %Loop through scales and angles to find best aligment
        for i_theta = 1:length(test_thetas)
            theta_i = test_thetas(i_theta);
            
            rot = [cos(theta_i) -sin(theta_i); sin(theta_i) cos(theta_i)];  
            vessel_xy_ori_i = vessel_xy_ori * exp(-2i*theta_i);

            for i_scale = 1:num_scales
                scale_i = args.scales(i_scale);

                %Transform vessel
                vessel_xy_i = bsxfun(@plus, vessel_xy*rot*scale_i, [offset_x offset_y]);

                %Sample orientations at this location
                patch_ori_i = interp2(ori_patch, vessel_xy_i(:,1), vessel_xy_i(:,2));
                
                %Compute difference between vessel and sampled orientations
                ori_diff = abs(angle(patch_ori_i .* vessel_xy_ori_i)/2);

                %Compute fit score - the sum of the sampled orientation magnitudes
                %(which we expect to be highest at the centre of the true
                %vessel) minus the sum of the orientation differences
                %(which should be minimised when the vessel is correctly
                %aligned)
                fit_score = -sum(ori_diff) + sum(abs(patch_ori_i)) ;

                %Check if this is the new best score
                if fit_score > best_fit
                    best_fit = fit_score;
                    best_theta = theta_i;
                    best_scale = scale_i;
                    best_xy = vessel_xy_i;
                    
                end
            end
        end
        vessel_xy = best_xy;
    end   

    if args.debug
        
        if i_ap <= 20
            figure;
            subplot(1,2,1); imgray(complex2rgb(ori_patch));
            plot(vessel_xy(:,1), vessel_xy(:,2), 'k', 'linewidth', 2);
            title(['Best \theta = ' num2str(round(180*best_theta/pi)) ', best scale =' num2str(best_scale), ', score =' num2str(best_fit)]);
            subplot(1,2,2); imgray(image_patch);
            plot(vessel_xy(:,1), vessel_xy(:,2), 'r', 'linewidth', 2);
            title(['Width at apex: ' num2str(vessel_width(ay, ax)) ', angle at apex: ' num2str(round(180*theta/pi))]);
            
        else
            break;
        end

    end
        
    %Save structure for this apex
    apex_candidate.vessel_xy = vessel_xy;
    apex_candidate.scale = scale*best_scale;
    apex_candidate.theta = theta+best_theta;
    apex_candidate.sr = sr;
    apex_candidate.sc = sc;
    apex_candidate.er = er;
    apex_candidate.ec = ec;
    apex_candidate.method = 'template_matching_aligned_g2d';
    apex_candidate.image_path = image_path;

    save([args.aam_dir 'apex' zerostr(i_ap, 4) '.mat'], 'apex_candidate');

    %Write out image patch
    imwrite(uint8(image_patch), [args.aam_dir 'images/candidate_apex' zerostr(i_ap,4) '.png']);

    %Write out a pts file we can read in to VXL
    fid1 = fopen([args.aam_dir 'points/candidate_apex' zerostr(i_ap,4) '.pts'], 'wt');
    fprintf(fid1, '%s \n', 'version: 1');
    fprintf(fid1, '%s %d \n', 'n_points:', num_pts);
    fprintf(fid1, '%s \n', '{'); 
    for i_p = 1:num_pts
        fprintf(fid1,'%.2f %.2f \n', vessel_xy(i_p,1), vessel_xy(i_p,2));
    end
    fprintf(fid1, '%s \n', '}');
    fprintf(fid1, 'nailfold: %s \n', image_path);
    fprintf(fid1, '%s %d \n', 'start_row:', sr);
    fprintf(fid1, '%s %d \n', 'start_col: ', sc);
    fclose(fid1);

    %Write entry for this image/pts pair in model .smd files
    str = ['candidate_apex' zerostr(i_ap,4) '.pts : candidate_apex' zerostr(i_ap,4) '.png'];
    fprintf(candidates_fid, '%s \n', str);
    
end

%Close up the candidates file
fprintf(candidates_fid, '%s \n', '}');
fclose(candidates_fid);
