function [failed_list] = convert_edge_markup_batch(contour_dir, selected_idx, display_output, overwrite)
%CONVERT_EDGE_MARKUP_BATCH *Insert a one line summary here*
%   [] = convert_edge_markup_batch(contour_dir)
%
% Inputs:
%      contour_dir - *Insert description of input variable here*
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
% Created: 16-Apr-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester

%Get list of input edge contour files
v_list = dir([contour_dir, '*vessel_edge.mat']);

if nargin < 2 || isempty(selected_idx)
    selected_idx = 1:length(v_list);
end
if nargin < 3
    display_output = 0;
end
if nargin < 4
    overwrite = 0;
end

%Loop through each file and convert the markup to compute the vessel
%contour
failed_list = [];
for i_ve = selected_idx
    
    if ~overwrite && exist([contour_dir v_list(i_ve).name(1:end-9) '_contour.mat'], 'file')
        display(['Vessel: ' v_list(i_ve).name ' contour already exists. Skipping (set overwrite to 1 if you want to recompute).']);
        continue;
    end
    
    %Load edge markup
    load([contour_dir v_list(i_ve).name]);
    
    %Process markup
    try
        [vessel_centre outer_edge inner_edge] = ...
            convert_edge_markup(outer_edge, inner_edge);
    catch err
        display(['Vessel: ' v_list(i_ve).name ' failed. Skipping.']);
        display(err.message);
        display('***');
        failed_list(end+1) = i_ve; %#ok
        continue;
    end
        
    
    %Save contour
    save([contour_dir v_list(i_ve).name(1:end-9) '_contour.mat'],...
        'vessel_centre', 'outer_edge', 'inner_edge');
    
    %Display output if requested
    if display_output
        figure; 
        axis ij equal; hold all;
        plot(...
            [outer_edge(:,1) inner_edge(:,1)]',...
            [outer_edge(:,2) inner_edge(:,2)]');

        plot(outer_edge(:,1),outer_edge(:,2), 'r-');
        plot(inner_edge(:,1),inner_edge(:,2), 'g-');    
        plot(vessel_centre(:,1),vessel_centre(:,2), 'k-', 'linewidth', 2);
    end
    
end
