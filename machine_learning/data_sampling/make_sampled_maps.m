function [] = make_sampled_maps(image_list, sampled_pts_list, save_dir, delete_pts_files)
%MAKE_SAMPLED_MAPS *Insert a one line summary here*
%   [] = make_sampled_maps(image_list, sampled_data, save_dir)
%
% Inputs:
%      image_list - *Insert description of input variable here*
%
%      sampled_pts_list - *Insert description of input variable here*
%
%      save_dir - *Insert description of input variable here*
%
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 30-Aug-2012
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester
if nargin < 4
    delete_pts_files = true;
end
num_images = length(image_list);
num_maps = length(sampled_pts_list);
d = floor(log10(num_images)+1);

for ii = 1:num_images
    
    %Load image to get its size
    im = u_load(image_list{ii,1});
    [rows cols dims] = size(im); %#ok
    clear im;
     
    sampled_maps = sparse( false(rows*cols, num_maps) );
    for jj = 1:num_maps
        sampled_pts = u_load(sampled_pts_list{jj});
        sampled_maps(sampled_pts{ii}, jj) = 1; %#ok
        
        if delete_pts_files && ii == num_images
            delete(sampled_pts_list{jj});
        end
    end
    image_list{ii,2} = [save_dir '/sampled_maps' zerostr(ii,d) '.mat'];
    save(image_list{ii,2}, 'sampled_maps');
end
save([save_dir '/sampled_maps_lookup.mat'], 'image_list');


    


