%%
clear;
clc;

filepath = '/home/mberks/isbe/datasets/MPIIGaze/Data/Normalized/p00/';
 
path = dir(filepath);
path = path(3:end);
files = {path.name};
 
total_num=0;
for num_f=1:length(files)
    readname = [filepath, files{num_f}];
    temp = load(readname);
    
    total_num = total_num+length(temp.filenames);
end

fprintf('File Path Ready!\n');
%%
Data=[];
Data.data = zeros(60,36,1, total_num*2);
Data.label = zeros(2, total_num*2);
Data.headpose = zeros(2, total_num*2);
Data.confidence = zeros(1, total_num*2);

index = 0;
 
for num_f=1:length(files)
    readname = [filepath, files{num_f}];
    temp = load(readname);
    
    num_data = length(temp.filenames(:,1));     
    for num_i=1:num_data
        % for left
        index = index+1;
        img = temp.data.left.image(num_i, :,:);
        img = reshape(img, 36,60);
        Data.data(:, :, 1, index) = img'; % filp the image
        
        Lable_left = temp.data.left.gaze(num_i, :)';
        theta = asin((-1)*Lable_left(2));
        phi = atan2((-1)*Lable_left(1), (-1)*Lable_left(3));
        Data.label(:,index) = [theta; phi];
 
        headpose = temp.data.left.pose(num_i, :);
        M = rodrigues(headpose);
        Zv = M(:,3);
        theta = asin(Zv(2));
        phi = atan2(Zv(1), Zv(3));
        Data.headpose(:,index) = [theta;phi];         
         
        % for right
        index = index+1;
        img = temp.data.right.image(num_i, :,:);
        img = reshape(img, 36,60);
        Data.data(:, :, 1, index) = double(flip(img, 2))'; % filp the image
         
        Lable_right = temp.data.right.gaze(num_i,:)';
        theta = asin((-1)*Lable_right(2));
        phi = atan2((-1)*Lable_right(1), (-1)*Lable_right(3));
        Data.label(:,index) = [theta; (-1)*phi];% flip the direction

        headpose = temp.data.right.pose(num_i, :); 
        M = rodrigues(headpose);
        Zv = M(:,3);
        theta = asin(Zv(2));
        phi = atan2(Zv(1), Zv(3));
        Data.headpose(:,index) = [theta; (-1)*phi]; % flip the direction
    end
    fprintf('%d / %d !\n', num_f, length(files)); 
end
 
fprintf('Saving\n');

Data.data = Data.data/255; %normalize
Data.data = single(Data.data); % must be single data, because caffe want float type
Data.label = single(Data.label);
Data.headpose = single(Data.headpose);

savename = 'MPII_data.h5';

store2hdf5(savename, Data.data, Data.label, 1); % the store2hdf5 function comes from https://github.com/BVLC/caffe/pull/1746
%% You can also use the matlab function for hdf5 saving:
%% hdf5write(savename,'/data', Data.data, '/label',[Data.label; Data.headpose]); 
fprintf('done\n');
