orig_list = dir('Z:\asymmetry_project\data\line_maps\2004_screening\abnormals\*.mat');
orig_names = get_mammo_info(orig_list);

%
temp_list = dir('Z:\asymmetry_project\data\linop_radial_maps2\2004_screening\abnormals\*016.mat');
temp_names = get_mammo_info(temp_list);
[set dummy idx] = intersect(temp_names, orig_names);
missing16 = setdiff(1:292, idx)';

temp_list = dir('Z:\asymmetry_project\data\linop_radial_maps2\2004_screening\abnormals\*032.mat');
temp_names = get_mammo_info(temp_list);
[set dummy idx] = intersect(temp_names, orig_names);
missing32 = setdiff(1:292, idx)';

temp_list = dir('Z:\asymmetry_project\data\linop_radial_maps2\2004_screening\abnormals\*064.mat');
temp_names = get_mammo_info(temp_list);
[set dummy idx] = intersect(temp_names, orig_names);
missing64 = setdiff(1:292, idx)';

temp_list = dir('Z:\asymmetry_project\data\linop_radial_maps2\2004_screening\abnormals\*128.mat');
temp_names = get_mammo_info(temp_list);
[set dummy idx] = intersect(temp_names, orig_names);
missing128 = setdiff(1:292, idx)';

temp_list = dir('Z:\asymmetry_project\data\linop_radial_maps2\2004_screening\abnormals\*256.mat');
temp_names = get_mammo_info(temp_list);
[set dummy idx] = intersect(temp_names, orig_names);
missing256 = setdiff(1:292, idx)';

rerun256 = setdiff(missing256, missing128);
rerun128 = setdiff(missing128, missing64);
rerun64 = setdiff(missing64, missing32);
rerun32 = setdiff(missing32, missing16);
rerun16 = missing16;

display(rerun256);
display(rerun128);
display(rerun64);
display(rerun32);
display(rerun16);





