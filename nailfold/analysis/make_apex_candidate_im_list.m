function [] = make_apex_candidate_im_list(image_names, label_dir, im_list_filename)
%F *Insert a one line summary here*
%   [] = f()
%
% Inputs:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 18-Jun-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

%Form full directory paths and create folder for HoGs
fid = fopen(im_list_filename, 'wt');

%Loop though each image
num_images = length(image_names);
for i_im = 1:num_images
    
    im_name = image_names{i_im} ;  
    display(['Processing image ' num2str(i_im) ', ' datestr(now)]);  
    
    load([label_dir im_name '_label.mat'], 'candidates_class');
    
    if isempty(candidates_class)
        continue;
    else
    
        num_candidates = size(candidates_class,1);

        for i_can = 1:num_candidates

            patch_name = [im_name '_' zerostr(i_can,3) '.png'];
            fprintf(fid, '%s %d\n', patch_name, candidates_class(i_can));
                    
        end
    end
end  
fclose(fid);
    
    
