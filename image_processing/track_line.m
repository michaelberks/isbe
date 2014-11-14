function [line_pts] = track_line(image_in, start_pt, varargin)
%TRACK_LINE *Insert a one line summary here*
%   [] = track_line(varargin)
%
% TRACK_LINE uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%
% Optional Arguments:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 24-Oct-2011
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin, '0', ...
    'first_dir', 1,...
    'dir_order', {'N', 'E', 'S', 'W', 'NE', 'SE', 'SW', 'NW'},...
    'max_length', inf);
clear varargin;

dir_order = circshift(args.dir_order, [0 args.first_dir-1]);
%init_dir = args.init_dir;

%Initialise line_pts with the given start point
line_pts = start_pt;

%Use while loop to keep tracking until no new points are found
go_on = true;
while go_on
   
   %Set go_on to false so we'll stop if we don't find a new point;
   go_on = false;
   
   %Check each direction in turn to find a new point
    for ii = 1:8
        dir = dir_order{ii};
       
        %Step in selected direction
        test_pt = move_step(line_pts(end,:), dir);

        %Check if point is on a line and not already visited
        if test_pt(1) > 0 && test_pt(1) <= size(image_in,1) && ...
            test_pt(2) > 0 && test_pt(2) <= size(image_in,2) && ...
                image_in(test_pt(1), test_pt(2)) && ~ismember(test_pt, line_pts, 'rows')
           
            %if so add it to the line_pts, set go_on to true and break out
            %the for loop
            line_pts = [line_pts; test_pt]; %#ok
           
            if size(line_pts, 1) < args.max_length
                go_on = true;
            end
            break;               
       end
   end
end


function [out_pt] = move_step(in_pt, step_dir)

    switch step_dir
        
        case 'N'
            out_pt = in_pt + [-1 0];
            
        case 'NE'
            out_pt = in_pt + [-1 1];
            
        case 'E'
            out_pt = in_pt + [0 1];
            
        case 'SE'
            out_pt = in_pt + [1 1];
            
        case 'S'
            out_pt = in_pt + [1 0];
            
        case 'SW'
            out_pt = in_pt + [1 -1];
            
        case 'W'
            out_pt = in_pt + [0 -1];
            
        case 'NW'
            out_pt = in_pt + [-1 -1];
    end
            