function synthesise_gaussian_wrapper(idx_in)

%set up arguments
% syn_args.PathToPyramidGM = [mberksroot, 'background/results/normal_2_single_gaussian'];
% syn_args.SamplePyramid = u_load([mberksroot, 'background/pyramid/normal_2/normal',...
%     zerostr(idx_in, 3), '.bmp_pyramid.mat']);
% syn_args.FilledImage = ones(size(syn_args.SamplePyramid{1,1}));
% syn_args.FilledImage(100:149, 100:149) = 0;
% syn_args.SaveFile = [mberksroot, 'background/syn/normal_2/single_gaussian',...
%     zerostr(idx_in, 3)];
% 
% %
% %profile synthesis
% profile on;
% mb_gmm_pyr_synthesis3(syn_args);
% profile off;
% profsave(profile('info'),[mberksroot, 'background/profiles/normal_2/single_gaussian',...
%     zerostr(idx_in, 3)]);

syn_args.SamplePyramid = u_load([mberksroot, 'background/pyramid/normal_2/normal006.bmp_pyramid.mat']);

r = 100;
[rows cols] = size(syn_args.SamplePyramid{1,1});
c_x = round(cols/2); c_y = round(rows/2);
[x y] = meshgrid(1:cols, 1:rows);
circ = (x - c_x).^2 + (y - c_y).^2 > r.^2;
clear r c_x c_y x y


syn_args.FilledImage = circ;
syn_args.ModelDir = [mberksroot, 'background/results/normal_2/'];
syn_args.ModelName = 'normal_2_model';
syn_args.CutOffLevel = 1;
syn_args.Plot = 0;

syn_args.SaveFile = [mberksroot, 'background/syn/normal006_normal_model'];
mb_gmm_pyr_synthesis(syn_args);
    