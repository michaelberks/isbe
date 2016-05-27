study_dir = 'C:\isbe\toyota\data\helen\';
anno_dir = [study_dir 'annotation\'];
annos_dir = [study_dir 'annotations\'];
image_dir = [study_dir 'images\'];
create_folder(annos_dir);
%%
for i_im = 1:2330
    filename = [anno_dir num2str(i_im) '.txt'];
    fid = fopen(filename,'r');
    frewind(fid);
    im_name = fgetl(fid);
    s = textscan(fid, '%.4f', 'commentstyle', '//', 'delimiter', ','); s = s{1};
    fclose(fid);
    fid = fopen([annos_dir im_name '.txt'], 'wt');
    fprintf(fid, '%.2f %.2f \n', s);
    fclose(fid);
end
%%
fid = fopen([study_dir 'trainnames.txt'],'r');
s = textscan(fid, '%s');
fclose(fid);
train_names = s{1};
fid = fopen([study_dir 'testnames.txt'],'r');
s = textscan(fid, '%s');
fclose(fid);
test_names = s{1};

%%
for i_im = 1:20
    im_name = train_names{i_im};
    extract_helen_patches(im_name, annos_dir, image_dir, [], [], 1);
end
%%
patch_dir = [study_dir 'image_patches\'];
create_folder(patch_dir);

fid = fopen([study_dir 'train_labels.txt'], 'wt');
for i_im = 1:length(train_names);
    im_name = train_names{i_im};
    extract_helen_patches(im_name, annos_dir, image_dir, patch_dir, fid, 0);
end
fclose(fid);

fid = fopen([study_dir 'test_labels.txt'], 'wt');
for i_im = 1:length(test_names);
    im_name = test_names{i_im};
    extract_helen_patches(im_name, annos_dir, image_dir, patch_dir, fid, 0);
end
fclose(fid);
%%
all_names = {train_names, test_names};
list_filenames = {[study_dir 'train_labels.txt'], [study_dir 'test_labels.txt']};
patch_dir = [study_dir 'image_eye_patches\'];
create_folder(patch_dir);

for i_type = 1:2
    fid = fopen(list_filenames{i_type}, 'wt');
    names_i = all_names{i_type};
    for i_im = 1:length(names_i);
        display(['Sampling from image ' num2str(i_im)]);
        im_name = names_i{i_im};
        extract_helen_eye_patches(im_name, ...
            'output_dir',           patch_dir,...
            'list_fid',             fid,...
            'n_samples_per_eye',    50,...
            'rand_sample_fun',      'uniform',...
            'rand_fun_params',      [32 32],...
            'patch_x',              -63.5:63.5,...
            'patch_y',              -63.5:63.5,...
            'do_left',              true,...
            'do_right',             true,...
            'mirror_right',         true,...
            'eye_l_idx',            135:154,...
            'eye_r_idx',            115:134,...
            'nose_idx',             42:58,...
            'mouth_idx',            59:87,...
            'iod_idx',              [135 115],...
            'base_iod',             100,...
            'plot',                 0);
    end
    fclose(fid);
end
