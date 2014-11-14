function [maxima_pos, maxima_vals, corr_map] = template_match_apexes(nailfold, templates, varargin)
%TEMPLATE_MATCH_APEXES *Insert a one line summary here*
%   [] = template_match_apexes(varargin)
%
% TEMPLATE_MATCH_APEXES uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%
% Optional Arguments:
%
% Outputs:
%       - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 06-Mar-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    'nailfold_mask', [],...
    'combine_method', 'product',...
    'sigma_1', 2,...
    'sigma_2', 2,...
    'feature_map',[],...
    'exclusion_zone', 20,...
    'threshold', 0, ...
    'plot', 0);

clear varargin;

%Make sure nailfold is not RGB 
nailfold = nailfold(:,:,1);

%Make a mask for the nailfold if one hasn't been supplied
if ~isempty(args.nailfold_mask)
    nailfold_mask = args.nailfold_mask;
    args = rmfield(args, 'nailfold_mask');
else
    nailfold_mask = make_nailfold_mosaic_mask(nailfold);
end

%Set up initial corr_map scores based on combination method
switch args.combine_method
    
    case 'product'       
        corr_map = ones(size(nailfold));
        combine_fun = @times;
    case 'sum'
        corr_map = zeros(size(nailfold));
        combine_fun = @plus;
        
    otherwise
        display(['Combine method ' args.combine method ' not recognised, using product']);
        corr_map = ones(size(nailfold));
        combine_fun = @times;
end

%Loop through eachtemplate
for i_te = 1:size(templates,1)
    
    switch templates{i_te,1}
        
        case 'intensity'
            feature_map = nailfold;
            
        case 'g1d'
            feature_map = gaussian_1st_derivative_gradient(nailfold, args.sigma_1);
            
        case 'g2d'
            feature_map = gaussian_2nd_derivative_line(nailfold, args.sigma_2);
            
        case 'vessel_prob'
            feature_map = args.feature_map; args = rmfield(args, 'feature_map');            
            
        otherwise
            display(['Template feature ' templates{i_te,1} ' not recognised, skipping this template']);
            continue;
            
    end

    %Apply normalised cross correlations using this template
    corr_map_i = mb_normxcorr2(templates{i_te,2}, feature_map);
    corr_map = feval(combine_fun, corr_map, corr_map_i);
    
    if args.plot
        figure; imgray(corr_map_i);
        figure; imgray(corr_map);
    end
end

%Search for local maxima
[maxima_pos, maxima_vals] = local_image_maxima(corr_map, args.exclusion_zone, nailfold_mask, args.threshold);
if args.plot   
    plot(maxima_pos(:,1), maxima_pos(:,2), 'rx');
end