function [] = initialise_aam_candidates_old(image_path, candidates_xy, varargin)

args = u_packargs(varargin, '0', ...
    'g2d_template', [],...
    'max_num_candidates', inf,...
    'thetas', -15:3:15,...
    'scales', 0.8:0.1:1.5,...
    'mean_shape', [],...
    'aam_dir', [],...
    'g2d_sigma', 4);
clear varargin;

%Load in nailfold
[~,~,file_ext] = fileparts(image_path);
if strcmp(file_ext, '.mat')
    nailfold = u_load(image_path);
else
    nailfold = imread(image_path); 
end
nailfold = nailfold(:,:,1);

%Compute the gaussian derivatives
[mag_2d] = gaussian_2nd_derivative_line(nailfold, args.g2d_sigma);

%--------------------------------------------------------------------------
%Now work out an initial alignment for each maxima and write out this
%candidate
num_candidates = min( size(candidates_xy,1), args.max_num_candidates);
candidate_scales = zeros(num_candidates,1);
candidate_thetas = zeros(num_candidates,1);

%Make a folder to store the candidates in
create_folder(args.aam_dir);
create_folder([args.aam_dir '\images\']);
create_folder([args.aam_dir '\points\']);
create_folder([args.aam_dir '\out_points\']);
candidates_fid = initialise_candidates_file(args.aam_dir);

%Precompute correlation denomiator for the template
x = repmat(-24:24, 49, 1);
y = repmat(-24:24, 49, 1)';
xy = [x(:) y(:)];

circle_mask = ~(x.^2 + y.^2 > 24^2);
N = sum(circle_mask(:));
template = args.g2d_template(circle_mask);

T = sum(template) / N;
T2 = sum(template.^2) / N;
denom_T = sqrt( max(T2 - T^2,0) );
num_pts = size(args.mean_shape, 1);

display(['Matching pose for ' num2str(num_candidates) ' apex candidates']);
for i_ap = 1:num_candidates
    display(['Matching pose for candidate ' num2str(i_ap)]);
    
    %Initial holder for maximum correlation score
    c_max = 0;
    
    %Loop through angles
    for theta = args.thetas
        rot = [cosd(theta) sind(theta); -sind(theta) cosd(theta)];

        %Loop through scales
        for scale = args.scales

            %Transform points given scale and angle and translate to
            %candidate position
            xya = xy * rot * scale;
            xa = reshape(xya(:,1) + candidates_xy(i_ap,1), 49, 49);
            ya = reshape(xya(:,2) + candidates_xy(i_ap,2), 49, 49);

            %Sample patch from 2nd deriv patch
            gd2_patch = interp2(mag_2d, xa, ya, 'bilinear');
            gd2_vec = gd2_patch(circle_mask);

            %Cross correlate with the template
            A = sum(gd2_vec) / N;
            A2 = sum(gd2_vec.^2) / N;

            denom = denom_T * sqrt( max(A2 - A^2,0) );
            numerator = (template' * gd2_vec) / N - A*T;

            c = numerator / denom;

            %Save this scale and angle if score the new maximum
            if c > c_max
                candidate_scales(i_ap) = scale;
                candidate_thetas(i_ap) = theta;
                vessel_xy = args.mean_shape * rot * scale;
                c_max = c;
            end
        end
    end

    %Compute position of RoI about candidate
    sr = max(1, floor(min(vessel_xy(:,2)) - 50 + candidates_xy(i_ap,2)));
    er = min(size(nailfold,1), ceil(max(vessel_xy(:,2)) + 50 + candidates_xy(i_ap,2)));
    sc = max(1, floor(min(vessel_xy(:,1)) - 50 + candidates_xy(i_ap,1)));
    ec = min(size(nailfold,2), floor(max(vessel_xy(:,1)) + 50 + candidates_xy(i_ap,1)));

    vessel_xy(:,1) = vessel_xy(:,1) - sc + candidates_xy(i_ap,1);
    vessel_xy(:,2) = vessel_xy(:,2) - sr + candidates_xy(i_ap,2);

    %Sample patch from image
    image_patch = nailfold(sr:er, sc:ec);

    %Save structure for this apex
    apex_candidate.vessel_xy = vessel_xy;
    apex_candidate.scale = candidate_scales(i_ap);
    apex_candidate.theta = candidate_thetas(i_ap);
    apex_candidate.sr = sr;
    apex_candidate.sc = sc;
    apex_candidate.er = er;
    apex_candidate.ec = ec;
    apex_candidate.method = 'template_matching_aligned_g2d';
    apex_candidate.image_path = image_path;

    save([args.aam_dir 'apex' zerostr(i_ap, 4) '.mat'], 'apex_candidate');
        
    %Write out image patch
    imwrite(uint8(image_patch), [args.aam_dir '/images/candidate_apex' zerostr(i_ap,4) '.png']);
    
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
