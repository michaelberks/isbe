bg_patch = imread('/home/mberks/matlab_data/background/images/mass/mass_bg_001.bmp');
filled_image = ones(size(bg_patch));
filled_image(200:299, 200:299) = 0;

profile on;
[synthesised_image] = mb_gmm_tex_synthesis('PathToTextureGMM',...
    '/home/mberks/matlab_data/background/models/mass_clustering_result',...
    'SeededImage', bg_patch, 'FilledImage', filled_image,...
    'SaveFrequency', 1000,...
    'SaveFile', '/home/mberks/matlab_data/background/syn/temp/syn_image'); 
profile off;
p = profile('info');
profsave(p,'/home/mberks/matlab_data/background/profiles/cluster_profile');
save('/home/mberks/matlab_data/background/syn/syn_bg_001', 'synthesised_image');
clear;
% path;
exit;
