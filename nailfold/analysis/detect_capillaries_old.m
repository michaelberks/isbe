function [] = detect_capillaries_old(image_path, varargin)
%DETECT_CAPILLARIES *Insert a one line summary here*
%   [] = detect_capillaries(varargin)
%
% DETECT_CAPILLARIES uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%
% Optional Arguments:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 14-Mar-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin, '0', ...
    'g2d_template', [],...
    'template_thresh', 0,...
    'max_num_candidates', 100,...
    'thetas', -15:3:15,...
    'scales', 0.8:0.1:1.5,...
    'mean_shape', [],...
    'aam_dir', [],...
    'aam_exe', 'ncm_sandpit_mb',...
    'aam_path', '',...
    'vessel_probs', [],...
    'g2d_sigma', 4, ...
    'g2d_vals', [],...
    'prob_thresh', 0.9,...
    'distal_angle', pi/4, ...
    'do_template_matching', 0,...
    'do_aam', 0, ...
    'do_final_vessel', 0, ...
    'delete_candidate_patches', 0,...
    'nailfold_mask', []);
clear varargin;

%Load in nailfold
nailfold = imread(image_path);   
nailfold = nailfold(:,:,1);

%Make a mask for the nailfold if one hasn't been supplied
if ~isempty(args.nailfold_mask)
    nailfold_mask = args.nailfold_mask;
    args = rmfield(args, 'nailfold_mask');
else
    nailfold_mask = make_nailfold_mosaic_mask(nailfold);
end

%Compute the gaussian derivatives
[mag_2d] = gaussian_2nd_derivative_line(nailfold, args.g2d_sigma);

if args.do_template_matching
%-------------------------------------------------------------------------
%Get potential apex candidates
templates = {'g2d', args.g2d_template};

[candidate_pos] = template_match_apexes(nailfold, templates,...
    'threshold', args.template_thresh, 'sigma_2', args.g2d_sigma);

display(['Detected ' num2str(size(candidate_pos,1)) ' apex candidates in initial template matching']);

%--------------------------------------------------------------------------
%Now work out an initial alignment for each maxima and write out this
%candidate
num_candidates = min( size(candidate_pos,1), args.max_num_candidates);
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
            xa = reshape(xya(:,1) + candidate_pos(i_ap,1), 49, 49);
            ya = reshape(xya(:,2) + candidate_pos(i_ap,2), 49, 49);

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
    sr = max(1, floor(min(vessel_xy(:,2)) - 50 + candidate_pos(i_ap,2)));
    er = min(size(nailfold,1), ceil(max(vessel_xy(:,2)) + 50 + candidate_pos(i_ap,2)));
    sc = max(1, floor(min(vessel_xy(:,1)) - 50 + candidate_pos(i_ap,1)));
    ec = min(size(nailfold,2), floor(max(vessel_xy(:,1)) + 50 + candidate_pos(i_ap,1)));

    vessel_xy(:,1) = vessel_xy(:,1) - sc + candidate_pos(i_ap,1);
    vessel_xy(:,2) = vessel_xy(:,2) - sr + candidate_pos(i_ap,2);

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
    imwrite(image_patch, [args.aam_dir '\images\candidate_apex' zerostr(i_ap,4) '.png']);
    
    %Write out a pts file we can read in to VXL
    fid1 = fopen([args.aam_dir 'points\candidate_apex' zerostr(i_ap,4) '.pts'], 'wt');
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

%--------------------------------------------------------------------------
end %do_template_matching
%--------------------------------------------------------------------------

if args.do_aam
%--------------------------------------------------------------------------
% Issue system command to fit AAM to each candidate
display('Fitting AAM');

candidates_path = [args.aam_dir 'all_candidates_test.smd'];
output_path = [args.aam_dir 'out_points\'];
quality_path = [args.aam_dir 'out_points\model_qualities.txt'];
cmd = [args.aam_exe ' -m ' args.aam_path ' -c ' candidates_path ' -o ' output_path ' -q ' quality_path];
system(cmd);

fid = fopen(quality_path);
q_txt = textscan(fid, '%s %f', 'delimiter', ':');
fclose(fid);
model_q = q_txt{2}; clear q_txt;
[sorted_model_qualities qidx] = sort(model_q, 'descend');


%--------------------------------------------------------------------------
% Now loop through each candidate again and add the model points and fit
% score to the apex structure
num_candidates = length(sorted_model_qualities);
for i_ap = 1:num_candidates
    j_ap = qidx(i_ap); 
    load([args.aam_dir 'apex' zerostr(j_ap, 4) '.mat'], 'apex_candidate');
    
    f1 = fopen([args.aam_dir 'out_points\candidate_apex' zerostr(j_ap,4) '.pts']);
    textscan(f1, '%[^{]');
    fgetl(f1);
    vessel_str = textscan(f1,'%[^}]');
    fclose(f1);
    apex_candidate.fitted_vessel_xy = str2num(vessel_str{1}{1}); %#ok
    apex_candidate.model_score = sorted_model_qualities(i_ap);
    
    %Save the new structure, now labelled in order of model score and
    %delete the old
    save([args.aam_dir 'apex' zerostr(i_ap, 4) '_aam.mat'], 'apex_candidate');
end

if args.delete_candidate_patches
    %Delete image, points and output points dir
    rmdir([args.aam_dir '\images\'], 's');
    rmdir([args.aam_dir '\points\'], 's');
    rmdir([args.aam_dir '\out_points\'], 's');
end

%--------------------------------------------------------------------------
end %do_aam
%--------------------------------------------------------------------------

if args.do_final_vessel
%--------------------------------------------------------------------------
% Finally, convert the 2nd deriv map into a CDF map of vessel probs
display('Fitting final vessels');

if isempty(args.vessel_probs)
    if isempty(args.g2d_vals)    
        args.g2d_vals = sort(mag_2d(nailfold_mask));
    end
    vessel_probs = interp1(args.g2d_vals, linspace(0,1,length(args.g2d_vals)), mag_2d);
    vessel_probs(~nailfold_mask) = 0;
else
    vessel_probs = args.vessel_probs;
    args = rmfield(args, 'vessel_probs');
end

num_candidates = length(dir([args.aam_dir 'apex*_aam.mat']));
for i_ap = 1:num_candidates
    
    %Load in apex structure
    load([args.aam_dir 'apex' zerostr(i_ap, 4) '_aam.mat'], 'apex_candidate');
    
    %Get patch from vessel_probs
    vessel_probs_patch = vessel_probs(...
        apex_candidate.sr:apex_candidate.er,...
        apex_candidate.sc:apex_candidate.ec);

    %Get vessels connected to the fitted apex
    apex_x = apex_candidate.fitted_vessel_xy(:,1);
    apex_y = apex_candidate.fitted_vessel_xy(:,2);
    vessel_i = bwselect(vessel_probs_patch > args.prob_thresh, apex_x , apex_y);
    
    if ~any(vessel_i(:))
        continue;
    end
    
    %Thin to get the centreline
    vessel_is = bwmorph(vessel_i, 'thin', inf);
    [skel_y skel_x] = find(vessel_is);
    
    %Use distance transform to compute width of vessel at each centre pixel
    width_map = bwdist(~vessel_i);
    skel_widths = width_map(vessel_is)*2 - 1;
    
    %Now project each fitted apex point to the centreline
    num_pts = length(apex_x);
    apex_widths = zeros(num_pts,1);
    apex_skel_x = zeros(num_pts,1);
    apex_skel_y = zeros(num_pts,1);
    for i_p = 1:num_pts
        [~, min_idx] = min((skel_x-apex_x(i_p)).^2 + (skel_y-apex_y(i_p)).^2);
        apex_widths(i_p) = skel_widths(min_idx);
        apex_skel_x(i_p) = skel_x(min_idx);
        apex_skel_y(i_p) = skel_y(min_idx);
    end
    apex_candidate.mean_apex_widths = mean(apex_widths);
    
%     %Finally, track down the skeleton from each end of the fitted apex
%     [left_pts] = track_line(vessel_is, [apex_skel_y(1) apex_skel_x(1)],...
%         'first_dir', 3, 'max_length', 100);
%     left_pts = [left_pts(end:-1:1,2) left_pts(end:-1:1,1)];
%     [right_pts] = track_line(vessel_is, [apex_skel_y(num_pts) apex_skel_x(num_pts)],...
%         'first_dir', 3, 'max_length', 100);
%     right_pts = [right_pts(:,2) right_pts(:,1)];
    
    %Final vessel 
%     apex_candidate.final_vessel = [left_pts; [apex_skel_x apex_skel_y]; right_pts];
    [vessel_y vessel_x] = find(vessel_i);
    apex_candidate.final_vessel = [vessel_x, vessel_y];
    
    centre_pt = (num_pts-1) / 2;
    global_x = apex_candidate.sc+apex_x(centre_pt);
    global_y = apex_candidate.sr+apex_y(centre_pt);
    if ~exist('distal_apexes', 'var');
        distal_apexes = [global_x global_y];
        apex_candidate.is_distal = true;
    else
        tan_vectors = abs(atan(...
            (distal_apexes(:,2)-global_y) ./...
            (distal_apexes(:,1)-global_x)));
        
        apex_candidate.is_distal = all(tan_vectors < args.distal_angle);
        if apex_candidate.is_distal
            distal_apexes = [distal_apexes; global_x global_y]; %#ok
        end
    end          
    save([args.aam_dir 'apex' zerostr(i_ap, 4) '_aam.mat'], 'apex_candidate');

end

end
%--------------------------------------------------------------------------
%End of main function
%--------------------------------------------------------------------------

