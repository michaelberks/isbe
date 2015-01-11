%Sorting out data for the 2 year study

%Notes:
% Base directory on MB's local machine is C:\isbe\nailfold\data\2_year_study
% Original data stored on the network drive at "N:\musculoskeletal\NCM only studies\2 yr data incl marina mark up".  
% Images (originally in .bmp) from dir "anonymous_bmap_files1" and its sub-dir "Additional mosaics 28-03-2014" 
% converted to .png/.mat format and stored in "anonymous_png"/"images" respectively. 
% .mat images are downsized by a factor of 2
function prepare_ncm_batch(local_root, varargin)

args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    'img_dir', 'images_bmp',...
    'img_format', 'bmp',...
    'b', []);
%% ------------------------------------------------------------------------
% 1) Convert and copy images from network directory, make a mosaic mask
%local_root = 'C:\isbe\nailfold\data\anniek\';

img_dir = [args.img_dir '\'];
mat_dir = 'images\';
mask_dir = 'fov_masks\';

create_folder([local_root mat_dir]);
create_folder([local_root mask_dir]);
%%
do_plot = 1;

im_list = dir([local_root img_dir '*.' args.img_format]);

num_ims = length(im_list);
for i_im = 1:num_ims

    display(['Preparing image ' num2str(i_im) ' of ' num2str(num_ims)]);
    im_name = im_list(i_im).name(1:end-4);

    %Load in image
    nailfold = imread([local_root img_dir im_name '.' args.img_format]);

    %Select single channel if in RGB format (R=G=B for grayscale bmps)
    if size(nailfold, 3) == 3
        nailfold = nailfold(:,:,1);
    end

    %Downsize and write as mat
    nailfold = imresize(nailfold, 0.5);
    save([local_root mat_dir im_name '.mat'], 'nailfold');

    %Make mosaic field-of-view mask
    fov_mask = make_nailfold_mosaic_mask(nailfold, 250, 5);
    save([local_root mask_dir im_name '_f_mask.mat'], 'fov_mask');

    %Show some images to check all is as expected
    if do_plot && i_im <= 5
        figure;
        subplot(2,1,1); imgray(nailfold);
        subplot(2,1,2); imgray(fov_mask);
    end
end


   