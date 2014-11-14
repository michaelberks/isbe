% GET_SPICULE_PARAMETERS - Collate spicule_data from set of masses into one structure 
%
% Usage: [mass_data spicule_data] = 
%       get_spicule_parameters(mass_files, s_width, n_pts, if_spic)
% 
% Arguments:
% mass_files - text file of newline separated file names to annotated mass
%  structures
% s_width - scalar integer specifying the half-width of the spicule normal
%  profile
% n_pts - scalar integer specifying the number of landmark points along the
%  length of the spicule
% if_spic - [0, 1], if set to 1, get spicule parameters, if 0 just masses
%
% Returns: 
% spicule_data - N-length structure, where N is the total number of spicules, 
% containing 5 fields:
%  w_vector - vector of parameters defining the spicule width
%  b_vector - vector of parameters defining the spicule brightness
%  o_vector - vector of offsets from the centre along the spicule length           
%  s_vector - m*2 matrix of (x y)-co-ordinates, where m is the number of 
%   landmark points, defining the spicule shape 
%  s_length - scalar of spicule length
%
% Notes:
% 
%
% See also: GET_SPICULE_DATA GENERATE_SPICULE_AM
%
% References:
%
% Author:   Michael Berks
%           Imaging Science and Biomedical Engineering
%           University of Manchester
%
function [mass_data spicule_data] = ...
    get_spicule_parameters(mass_files, s_width, n_pts, if_spic)

    if nargin < 4
        if_spic = 1;
    end
    
    fid = fopen(mass_files);
    file_names = textscan(fid, '%s');
    clear fid mass_files;
    
    ns = 1; %running count of the number of spicules
    
    for ii = 1:size(file_names{1}, 1)
        
        [sd md] = get_spicule_data(file_names{1}{ii}, s_width, n_pts);
        mass_data(ii) = md;
        if if_spic
            for jj = 1:length(sd)
                %
                % Extract width, offset and brightness vectors from sd
                %
                pd = fit_ellipse_profile(sd(jj).profile, 30, 10, 4, n_pts);
                spicule_data(ns).w_vector = pd.w_params;
                spicule_data(ns).b_vector = pd.b_params;
                spicule_data(ns).o_vector = pd.x0_local;
                spicule_data(ns).s_vector = sd(jj).s_shape;
                spicule_data(ns).s_length = sd(jj).s_length;
                spicule_data(ns).s_orientation = sd(jj).s_orientation;
                spicule_data(ns).s_location = sd(jj).s_location;
                spicule_data(ns).s_distance = sd(jj).s_distance;
                ns = ns + 1;
            end
        end
    end
    if ns == 1
        spicule_data = struct;
    end
end