mb_build_pyr_tsvq('PyramidDir', [mberksroot, 'background/pyramid/normal_2/'],...
    'OutputDir', [mberksroot, 'background/trees/normal_2/']);

% clear;
% 
% pyramid1 = u_load('C:\isbe\dev\background\pyramid\normal_2\normal001.bmp_pyramid.mat');
% %pyramid2 = u_load('C:\isbe\dev\background\pyramid\normal_2\normal002.bmp_pyramid.mat');
% filled_image = logical(ones(size(pyramid1{2,1})));
% filled_image(100:249, 100:299) = 0;
% wei_args.FilledImage = filled_image;
% wei_args.SampleData = u_load('C:\isbe\dev\background\trees\normal_2\normal002_tree\tsvq_names.mat');
% wei_args.TargetPyramid = pyramid1;
% wei_args.WindowSize1 = 11;
% wei_args.WindowSize2 = 7;
% 
% %%
% % wei_args.K_method = 'standard';
% % profile on
% % [syn_im_full, new_pyr_full] =  mb_wei_pyr_synthesis(wei_args);
% % profile off
% % profsave(profile('info'),'C:\isbe\dev\background\profiles\wei_full')
% 
% wei_args.K_method = 'tree';
% profile on
% [syn_im_tree, new_pyr_tree] =  mb_wei_pyr_synthesis(wei_args);
% profile off
% profsave(profile('info'),'C:\isbe\dev\background\profiles\wei_tree');
% figure; imagesc(syn_im_tree); colormap(gray(256)); axis image;
% 
% %% 
% clear
% pyramid2 = u_load([mberksroot, 'background/pyramid/normal_2/normal002.bmp_pyramid.mat']);
% 
% dirname = [mberksroot, 'background/trees/normal_2/normal002_tree/'];
% mkdir(dirname);
% 
% %profile on;
% 
% tsvq_names = cell(7, 5, 2);
% for level = 2:5
%     for ori = 1:5
% %         k_args.Image1 = pyramid2{level+1, ori};
% %         k_args.Image2 = pyramid2{level, ori};
% %         k_args.WindowSize1 = 11;
% %         k_args.WindowSize2 = 7;
% %         k_args.Method = 'tree';
% %         k_args.BothLevels = false;
% % 
% %         tsvq_tree = mb_knearest_build(k_args); %#ok
% 
%         tree_name = [dirname, 'tsvq_', num2str(level), '_', num2str(ori), '_1'];
% %         save(tree_name, 'tsvq_tree');
%         tsvq_names{level,ori, 1} = tree_name;
% %         clear tsvq_tree;
% %         
% %         k_args.BothLevels = true;
% %         tsvq_tree = mb_knearest_build(k_args);
%         tree_name = [dirname, 'tsvq_', num2str(level), '_', num2str(ori), '_2'];
% %         save(tree_name, 'tsvq_tree');
%         tsvq_names{level,ori, 2} = tree_name;
% %         clear k_args tsvq_tree;
%     end
% end
% 
% save([dirname, 'tsvq_names'], 'tsvq_names');
% % profile off;
% % profsave(profile('info'),[dirname, '/wei_profile']);
