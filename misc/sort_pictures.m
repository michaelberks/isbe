function sort_pictures(pic_dir, pic_name)

if ~strcmp(pic_dir(end), filesep)
    pic_dir = [pic_dir filesep];
end

%get list of JPG files in directory
pic_list = dir([pic_dir '*.jpg']);

num_pics = length(pic_list);
num_length = length(num2str(num_pics));

%get dates of each picture
for ii = 1:num_pics
    pic_info = imfinfo([pic_dir pic_list(ii).name]);
    pic_dates(ii,:) = pic_info.DateTime; %#ok
end

%sort picture
[dummy,idx] = sortrows(pic_dates);
%rename pictures
for ii = 1:num_pics
    movefile(...
        [pic_dir pic_list(idx(ii)).name],...
        [pic_dir pic_name zerostr(ii,num_length) '.jpg']);
end