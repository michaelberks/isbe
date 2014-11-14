function [] = fit_aam_to_candidates_temp(image_path, candidates_xy, varargin)

args = u_packargs(varargin, '0', ...
    'g2d_template', [],...
    'max_num_candidates', inf,...
    'thetas', -15:3:15,...
    'scales', 0.8:0.1:1.5,...
    'mean_shape', [],...
    'aam_dir', [],...
    'aam_exe', 'ncm_sandpit_mb',...
    'aam_path', '',...
    'g2d_sigma', 4, ...
    'delete_candidate_patches', 0);
clear varargin;

%Load in nailfold
[~,~,file_ext] = fileparts(image_path);
if strcmp(file_ext, '.mat')
    nailfold = u_load(image_path);
else
    nailfold = imread(image_path); 
end
nailfold = nailfold(:,:,1);

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

display(['Matching pose for ' num2str(num_candidates) ' apex candidates']);
for i_ap = 1:num_candidates
    apex_candidate = u_load([args.aam_dir 'apex' zerostr(i_ap, 4) '.mat']);

    %Sample patch from image
    image_patch = nailfold(apex_candidate.sr:apex_candidate.er, apex_candidate.sc:apex_candidate.ec);   
    vessel_xy = apex_candidate.vessel_xy;
      
    %Write out image patch
    imwrite(uint8(image_patch), [args.aam_dir '\images\candidate_apex' zerostr(i_ap,4) '.png']);
    
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
    fprintf(fid1, '%s %d \n', 'start_row:', apex_candidate.sr);
    fprintf(fid1, '%s %d \n', 'start_col: ', apex_candidate.sc);
    fclose(fid1);
    
    %Write entry for this image/pts pair in model .smd files
    str = ['candidate_apex' zerostr(i_ap,4) '.pts : candidate_apex' zerostr(i_ap,4) '.png'];
    fprintf(candidates_fid, '%s \n', str);
end

%Close up the candidates file
fprintf(candidates_fid, '%s \n', '}');
fclose(candidates_fid);

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
    apex_candidate.fitted_vessel_xy = str2num(vessel_str{1}{1});
    apex_candidate.model_score = sorted_model_qualities(i_ap);
    
    %Save the new structure, now labelled in order of model score and
    %delete the old
    save([args.aam_dir 'apex' zerostr(i_ap, 4) '_aam.mat'], 'apex_candidate');
end

if args.delete_candidate_patches
    %Delete image, points and output points dir
    movefile(quality_path, [args.aam_dir 'model_qualities.txt']);
    rmdir([args.aam_dir '\images\'], 's');
    rmdir([args.aam_dir '\points\'], 's');
    rmdir([args.aam_dir '\out_points\'], 's');
end