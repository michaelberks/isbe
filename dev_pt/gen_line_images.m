addpath(genpath('S:\projects\mammography\matlab\'));

generate_synthetic_curve_images(...
	'num_images',100,...
	'save_dir','',... % the mandatory arguments
	'bg_dir', 'U:\projects\mammography\data\synthetic_backgrounds\smooth512\train\',...
	'bg_fmt', 'mat',...
	'save_path', 'U:\projects\mammography\data\synthetic_lines2\',...
	'save_images',true,...
  'contrast_range', [1 8],...
	'plot',0);
