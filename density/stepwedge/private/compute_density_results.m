function [density_results] = ...
            compute_density_results(gland_thickness, breast_thickness,...
                taper_region_idx, fid1, fid2, fid3, fid4, filename)
      
%DENSITY_RESULTS *Insert a one line summary here*
%   [] = density_results(gland_thickness,breast_thickness,taper_region_idx)
%
% Inputs:
%      gland_thickness- *Insert description of input variable here*
%
%      breast_thickness- *Insert description of input variable here*
%
%      taper_region_idx- *Insert description of input variable here*
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
% Created: 12-May-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

% set area of pixel value (in cm^2)
% image has not been resized and has a 250um pixel size
pixel_area = 0.025*0.025;

% work out area and volume of breast based on breast mask and thickness

% work out area and volume of breast and glandular tissue
area_b = sum(breast_thickness(:) > 0)*pixel_area; % this is in units of cm^2
area_g = sum(gland_thickness(:) > 0)*pixel_area;

volume_b = sum(breast_thickness(:))*pixel_area*0.1; % units are cm^3 
volume_g = sum(gland_thickness(:))*pixel_area*0.1; 

dense_by_area_1 = 100 * area_g / area_b ;
dense_by_vol_1  = 100 * volume_g / volume_b ;

% Area and volume %density will be VERY different because the area
% measure of gland includes all pixels where gland_thickness>0
% There will be very few pixels where composition is 100% fat!
% Therefore, include a number of different measures:
% area_g5, area_g10... where the number denotes gland thickness as a
% proportion of breast thickness

area_g5 = sum(gland_thickness(:)./breast_thickness(:) > 0.05)*pixel_area;
area_g10 = sum(gland_thickness(:)./breast_thickness(:) > 0.10)*pixel_area;
area_g15 = sum(gland_thickness(:)./breast_thickness(:) > 0.15)*pixel_area;
area_g20 = sum(gland_thickness(:)./breast_thickness(:) > 0.20)*pixel_area;
area_g25 = sum(gland_thickness(:)./breast_thickness(:) > 0.25)*pixel_area;

density_area_g5 = 100 * area_g5 / area_b;
density_area_g10 = 100 * area_g10 / area_b;
density_area_g15 = 100 * area_g15 / area_b;
density_area_g20 = 100 * area_g20 / area_b;
density_area_g25 = 100 * area_g25 / area_b;

density_results.volume_b = volume_b;
density_results.volume_g = volume_g;
density_results.area_b = area_b;
density_results.dense_by_area_1 = dense_by_area_1; 
density_results.dense_by_vol_1 = dense_by_vol_1;
density_results.density_area_g5 = density_area_g5;
density_results.density_area_g10 = density_area_g10;
density_results.density_area_g15 = density_area_g15;
density_results.density_area_g20 = density_area_g20;
density_results.density_area_g25 = density_area_g25;
density_results.area_g = area_g;
density_results.area_g5 = area_g5; 
density_results.area_g10 = area_g10;
density_results.area_g15 = area_g15;
density_results.area_g20 = area_g20;
density_results.area_g25 = area_g25;
density_results.max_thickness = max(breast_thickness(:));
density_results.min_thickness = min(breast_thickness(:));

%If we are given file IDs, then print the information to the text files
if nargin > 3
    fprintf(fid1,'%s   %7.4f   %7.4f   %7.4f \n',filename, area_b, area_g, dense_by_area_1);
    fprintf(fid2,'%s   %7.4f   %7.4f   %7.4f \n',filename, volume_b, volume_g, dense_by_vol_1);
    fprintf(fid3,'%s   %7.4f   %7.4f   %7.4f   %7.4f   %7.4f   %7.4f   %7.4f   %7.4f    %7.4f     %7.4f     %7.4f     %7.4f \n',...
        filename, area_g, area_g5, area_g10, area_g15, area_g20, area_g25, dense_by_area_1,...
        density_area_g5, density_area_g10, density_area_g15, density_area_g20, density_area_g25 );
    fprintf(fid4,'%s   %7.4f   %7.4f \n', filename, max(breast_thickness(:)), min(breast_thickness(:)) );
end
