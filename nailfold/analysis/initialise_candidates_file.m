%--------------------------------------------------------------------------
function fid = initialise_candidates_file(aam_dir)

    fid(1) = fopen([aam_dir 'all_candidates_test.smd'], 'wt');
    
    fprintf(fid, '%s \n', '// Text file describing test data for vessel apex models');
    fprintf(fid, '\n');
    fprintf(fid, '%s \n', '// Directory containing images');
    fprintf(fid, '%s \n', ['image_dir: ' aam_dir 'images\']);
    fprintf(fid, '%s \n', '// Directory containing points');
    fprintf(fid, '%s \n', ['points_dir: ' aam_dir 'points\']);
    fprintf(fid, '\n');
    fprintf(fid, '%s \n', '// Details of points : images');
    fprintf(fid, '\n');
    fprintf(fid, '%s \n', 'training_set:');
    fprintf(fid, '%s \n', '{');