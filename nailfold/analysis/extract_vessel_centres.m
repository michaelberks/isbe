function [vessel_centre vessel_prob vessel_ori vessel_width] = ...
    extract_vessel_centres(vessel_prob, vessel_ori, vessel_width, varargin)
%EXTRACT_VESSEL_CENTRES *Insert a one line summary here*
%   [] = extract_vessel_centres()
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
args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    'prob_sigma',           2,...
    'ori_sigma',            2,...
    'width_sigma',            2,...
    'strong_vessel_thresh', 0.25,...
    'weak_vessel_thresh',   0);
    
%Smooth the vessel probs
if args.prob_sigma
    g_prob = gaussian_filters_1d(args.prob_sigma);
    g_prob = g_prob / sum(g_prob);    
    vessel_prob = conv2(g_prob', g_prob, vessel_prob, 'same');
end
if args.ori_sigma
    g_ori = gaussian_filters_1d(args.ori_sigma);
    g_ori = g_ori / sum(g_ori);
    vessel_ori = conv2(g_ori', g_ori, vessel_ori, 'same');
end 
if args.width_sigma
    g_width = gaussian_filters_1d(args.width_sigma);
    g_width = g_width / sum(g_width);
    vessel_width = conv2(g_width', g_width, vessel_width, 'same');
end 

%Compute NMS centrelines
vessel_nms = mb_non_maximal_supp(vessel_prob, angle(vessel_ori)/2);
strong_vessels = vessel_nms > args.strong_vessel_thresh;
if any(strong_vessels(:))
    [rstrong cstrong] = find(strong_vessels);
    vessel_centre_mask = bwselect(vessel_nms > args.weak_vessel_thresh, cstrong, rstrong, 8);
else
    vessel_centre_mask = strong_vessels;
end
[vessel_centre.y vessel_centre.x] = find(vessel_centre_mask); 
vessel_centre.prob = vessel_prob(vessel_centre_mask);
vessel_centre.ori = vessel_ori(vessel_centre_mask); 
vessel_centre.width = vessel_width(vessel_centre_mask);    
    
    
