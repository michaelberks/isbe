% clc; 
clear;

% create ideal orientation map
imsz = 100; mass_c2 = (imsz+1)*ones(1,2);
[xx,yy] = meshgrid(-imsz:imsz,-imsz:imsz);
tt0 = atan2(-yy,xx);

% choose noise levels
sig_list = logspace(0,-2,10);
sig_list = 0.0;
ntrials = 1;

stats = zeros(ntrials,length(sig_list));
for isig = 1:length(sig_list)
	for trial = 1:ntrials
		% add noise to ideal orimap
		sig = sig_list(isig);
		orimap = mod( tt0+sig*randn(size(tt0)),pi );
		linemap = (rand(size(orimap))>0.0); % change this threshold if need be

		map_sz = 0; spacing = 4;
% 		r_max = 80; r_min = r_max/2; R = r_min/4;
		r_max = 50; r_min = r_max/2; R = r_min/2;
		[f1,f2,k_mask] = karssemeijer_radial_projection_multiscale(...
							linemap,orimap,...
							'num_angles',1,...
							'r_min',r_min,'r_max',r_max,'R',R,...
							'x_min',mass_c2(1)-map_sz,...
							'x_max',mass_c2(1)+map_sz,...
							'y_min',mass_c2(2)-map_sz,...
							'y_max',mass_c2(2)+map_sz,...
							'spacing',spacing);

		stats(trial,isig) = f1;
		
% 		k_map_f1	= nan(size(k_mask));
% 		k_map_f1(k_mask) = f1(:,1);
% 		k_map_f1	= k_map_f1( mass_c2(2)-map_sz:spacing:mass_c2(2)+map_sz,...
% 								mass_c2(1)-map_sz:spacing:mass_c2(1)+map_sz );
% 		k_map_f2	= nan(size(k_mask));
% 		k_map_f2(k_mask) = f2(:,1);
% 		k_map_f2	= k_map_f2( mass_c2(2)-map_sz:spacing:mass_c2(2)+map_sz,...
% 								mass_c2(1)-map_sz:spacing:mass_c2(1)+map_sz );
% 
% 		cmap	= [0 0 0; hsv(255)];
% 		figure(1); clf; colormap(cmap); sbsz = [2,3];
% 			image(uint8(1+normim(tt)*255));
% 			axis('image');
% 		figure(2); clf;
% 		subplot(1,2,1); imagesc(k_map_f1); axis('image');
% 		subplot(1,2,2); imagesc(k_map_f2); axis('image');
	end
end

figure(10);
% 	stats1 = reshape(stats(1,:,:),[length(sig_list),ntrials])';
	stats1 = stats;
	boxplot(stats1);
	set(gca,'xticklabel',num2str(sig_list','%5.2f'));
	xlabel('Noise Sigma (rad)');
	
figpath = 'U:\projects\mammography\figs\kmaps\';
set(gcf,'paperposition',[0,0,16,10]);
% exportfig([figpath,'karssemeijer1']);

