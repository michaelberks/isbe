function [] = fit_aam_to_candidates(varargin)

args = u_packargs(varargin, '0', ...
    'aam_dir', [],...
    'aam_exe', 'ncm_sandpit_mb',...
    'aam_path', '',...
    'delete_candidate_patches', 0,...
    'sort_by_aam_fit', 0);
clear varargin;

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

num_candidates = length(model_q);
if args.sort_by_aam_fit
    [sorted_model_qualities qidx] = sort(model_q, 'descend');
else
    qidx = 1:num_candidates;
    sorted_model_qualities = model_q;
end

%--------------------------------------------------------------------------
% Now loop through each candidate again and add the model points and fit
% score to the apex structure

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
    if args.delete_candidate_patches
        %Delete image, points and output points dir
        delete([args.aam_dir 'apex' zerostr(j_ap, 4) '.mat']);
    end
end

if args.delete_candidate_patches
    %Delete image, points and output points dir
    movefile(quality_path, [args.aam_dir 'model_qualities.txt']);
    rmdir([args.aam_dir '\images\'], 's');
    rmdir([args.aam_dir '\points\'], 's');
    rmdir([args.aam_dir '\out_points\'], 's');
end