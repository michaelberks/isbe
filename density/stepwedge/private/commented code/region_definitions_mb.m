function [ region_lims ] = region_definitions_mb(im_sz, region_selection, large_image, left_breast, padding)
%REGION_DEFINITIONS : Depending on the size of the image and which set of regions required; outputs the marker regions. Each region is defined by
%					  a number to denote which specific region it is. For example a small image will have 6 regions, 3 on each side. Each of the 3
%					  regions on a side will be an area of the number 1, 2 or 3 with the corresponding area on the other side of the image
%					  being a region of the same number. This is done so that markers found in a region can be identified as being found in a 
% 					  specific region so that the discrimination is made easier. Any non-region is set to zero.
% 
%
% Inputs:
%			image				image that you want to find markers in [2D array]
%			region_selection	defines which regions to output for small images [#(1 or 2)]
%
% Outputs:
%			region_map			Map of the marker regions in the image. Marker regions are areas of one's with the rest of the image being zeros
%								[2D array]
%
% Example:
%
% Notes:	Large image regions already have an added 50 pixel padding on all edges and so are written differently to the small regions in which 
%			the padding is added in this script. Hence the difference in the code for the two. 
%
%			
% See also:
%
% Created: 18-06-2013
% Author: Euan Allen
% Email : euan.allen@gmail.com

%Set default padding
if ~exist('padding', 'var')
    if large_image
        padding = 0;
    else
        padding = 50;
    end
end

%Get size of image
rows = im_sz(1);
cols = im_sz(2);

if large_image
    region_lims = [ 1037 1681,    4  612;
                    2090 2644,    4  612;
                    3105 3620,    4  612;
                    4116 4620,    4  612;
                    1037 1681, 6175 6682;
                    2090 2644, 6175 6682;
                    3105 3620, 6175 6682;
                    4116 4620, 6175 6682];    
else
    if region_selection == 1
    
%         region_lims = [    1  470,  920 1314;
%                            1  453, 1977 2263;
%                            1  423, 2978 3286;                       
%                         4909 5202,  943 1362;
%                         4920 5231, 1941 2314;
%                         4969 5267, 2938 3277];

        region_lims = [  350  750,       51  450;
                        1350 1650,       51  450;
                        2350 2650,       51  450;
                         350  750, rows-399 rows;
                        1350 1650, rows-399 rows;
                        2350 2650, rows-399 rows];
                    

    else %if region_selection ==2
    
%         region_lims = [  149  469,  373  732;
%                           65  468, 1345 1676;
%                            1  466, 2364 2653;                       
%                         4917 5254,  409  739;
%                         4924 5307, 1263 1669;
%                         4932 5353, 2310 2659];
        region_lims = [  900 1300,       1  400;
                        2000 2300,       1  400;
                        3000 3300,       1  400;
                         900 1300, rows-399 rows;
                        2000 2300, rows-399 rows;
                        3000 3300, rows-399 rows];
    end
    
end

if left_breast
    region_lims(:,[1 2]) = cols - region_lims(:,[2 1]) + 1;
end
    
%Apply any padding
region_lims(:,[1 3]) = region_lims(:,[1 3]) - padding;
region_lims(:,[2 4]) = region_lims(:,[2 4]) + padding;

%Make sure the region limits lie within the boundaries of this specific
%image
region_lims(region_lims(:,1)<1,1) = 1;
region_lims(region_lims(:,3)<1,3) = 1;
region_lims(region_lims(:,2)>cols,2) = cols;
region_lims(region_lims(:,4)>rows,4) = rows;