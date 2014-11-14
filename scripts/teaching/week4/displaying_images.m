%--------------------------------------------------------------------------
% ------------------------- Displaying Images -----------------------------
%--------------------------------------------------------------------------

%Script showing examples of how the main imaging functions in matlab work

%--------------------------------------------------------------------------
% 2D arrays using colormaps
%--------------------------------------------------------------------------
data_dir = 'P:\MatlabTutorial\data\';
%load in an image
load([data_dir '024RCC.mat']); %loads in a data structure mammogram

%mammo is uint8, lets also make a double copy of it
mammogram_d = double(mammogram);

%-----------------------------------------
% Using image
%-----------------------------------------
%%
%Display it using the most basic imaging function 'image'
figure; 
image(mammogram);
%Image uses colormap 'jet' and is 'squashed' to fit in the screen
%%
figure; 
image(mammogram);
axis image;
colormap gray;
%%
figure; 
image(mammogram_d*2);
axis image;
colormap gray;
%We've doubled all the values in the image - many now lie above the max
%color range (256) and appear as white
%%
figure; 
image(mammogram_d/2);
axis image;
colormap gray;
%We've halved all the values in the image, the display is now twice as dark
%%
figure; 
image(-mammogram_d);
axis image;
colormap gray;
%By taking the negative of the image, all the values now lie outside the
%colour range, the image is black
%%
%------------------------------------
%Using imagesc
%------------------------------------
%Repeat the above images, but using imagesc, which auto-scales the image
%colour ranges
%Note we still need to use the axis image and colormap gray functions...
figure; imagesc(mammogram);
figure; imagesc(mammogram); axis image; colormap gray;
figure; imagesc(mammogram_d*2); axis image; colormap gray;
figure; imagesc(mammogram_d/2); axis image; colormap gray;
figure; imagesc(-mammogram_d); axis image; colormap gray;
%%
%------------------------------------
%Using imshow
%------------------------------------
figure; imshow(mammogram);
figure; imshow(mammogram_d);
figure; imshow(mammogram_d/2);
figure; imshow(-mammogram_d);
%%
%--------------------------------------------------------------------------
% RGB arrays
%--------------------------------------------------------------------------
%Now see what happen with RGB arrays
cancer_map = imread([data_dir 'f_map.jpg']);

%Cancer map is read is an RGB array of uint8s, again, take a double copy
cancer_map_d = double(cancer_map);

%-----------------------------------------
% Using image
%-----------------------------------------
%%
figure; image(cancer_map); axis image; %Image is displayed in color
figure; image(cancer_map); axis image; colormap gray; %Image is displayed in color - changing the colormap makes no difference
%%
figure; image(cancer_map_d); axis image; 
%Causes error, because image is a double, it's values must lie in the range
%0, 1
%%
%Dividing by 255 fixes this
figure; image(cancer_map_d/255); axis image;

%Dividing in half again makes the image twice as dark
figure; image(0.5*cancer_map_d/255); axis image;
%%
%-----------------------------------------
% Using imagesc - makes no difference, we're not scaling a colormap, we're
% interpreting raw values in the input array
%-----------------------------------------
figure; imagesc(cancer_map); axis image;
figure; imagesc(cancer_map_d/255); axis image;
figure; imagesc(0.5*cancer_map_d/255); axis image;
figure; imagesc(cancer_map_d); axis image; %Still causes error
%%
%-----------------------------------------
% Using imshow
%-----------------------------------------
figure; imshow(cancer_map);
figure; imshow(cancer_map_d/255);
figure; imshow(0.5*cancer_map_d/255);
figure; imshow(cancer_map_d); %No longer causes an error, but outside the range are scaled to 0 or 1
