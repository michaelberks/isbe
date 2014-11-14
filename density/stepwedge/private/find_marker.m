function [xc, yc, rad] = find_marker(mammo, xi, yi, reject, marker_fig)
%FIND_MARKER given an initial approximation to a breastt hickness marker
%centre, find the precise centre automatically and estimate the radius of
%the marker
%
% Arguments:
%
% 'mammo'
%   - Mammogram containing breast thickness markers
%
% 'xi,yi'
%   - Initial estimate of marker centre
%
% 'reject'
%   - if 0 (i.e. for first attempt) marker will be found fully
%   automatically. If 1 (i.e. user has marker first attempt as reject),
%   user is asked to select points on edge to detect marker
%
% 'marker_fig'
%   - handle to the figure of marker images (e.g. generated in calling function
%   marker detect) so we can draw ouput in this figure
%
%
% Outputs:
%
% 'xc, yc'
%   - x,y coordinates of marker centre
%
% 'rad'
%   - radius of marker
%
% See also: STEPWEDGE, MARKERDETECT
%
% Created: 09-May-2006
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 5125 
% Copyright: (C) University of Manchester
    
    %get the dimensions of the image
    [r c] = size(mammo);
    
    %specify size of box in which to search for marker
    box_size = 100;
    
    % get boundaries for box, ensuring we don't sample outside the image      
    x_start = max(1, xi-box_size);
    x_end = min(c, xi+box_size);
    y_start = max(1, yi-box_size);
    y_end = min(r, yi+box_size);

    %adjust the co-ordinates of the original point
    xi = xi - x_start + 1;
    yi = yi - y_start + 1;

    % extract the region around this point (processing would be too
    % slow to operate on whole image)
    small = mammo(y_start:y_end, x_start:x_end);
    
    %Set maximum radius of valid centres
    max_rc = 40;
        
    %Unless this is a rejection, try and find the marker automatically
    if ~reject
   
        [rs cs] = size(small);
        
        %median filter the region to remove noise
        small_med = medfilt2(small, [7 7]);

        %Work out valid centres
        xx = repmat(1:cs, rs, 1);
        yy = repmat((1:rs)', 1, cs);
        cc = (xx-xi).^2 + (yy-yi).^2 < max_rc^2;
        cc(1:2:end, 1:2:end) = 0;
        [centres_y centres_x] = find(cc);

        %Do an initial canny edge detector to work out strong straight lines so
        %we can discard these
        cannysmall = canny_edge(small_med,0.5,2,[],0.5);

        %Compute the initial ignore_map, discarding outside the large central
        %circle
        min_rad = 20;
        max_rad = 75;
        ignore_map = (xx-xi).^2 + (yy-yi).^2 > max_rad.^2;

        %separate edges into 8-connected structures
        [edge_labels n] = bwlabel(cannysmall,8);

        %for each edge, add straight line segments to the ignore map
        for kk = 1:n
            these_edges = edge_labels == kk;

            %if it's of reasonable size
            if sum(these_edges(:)) > 20

                %work out the gradient
                [ey ex] = find(these_edges);
                [dummy fxy] = gradient([ex ey]);
                std_fxy = std(fxy);

                %if this is fairly constant (i.e. has a low SD) then we
                %probably have a straight line
                if sum(std_fxy) < 1;
                    ignore_map(these_edges) = 1;
                end
            end
        end

        %Now repeat the cann edge detector using the new ignore map
        circular_edges = canny_edge(small_med,[],2,0.995,0.5, ignore_map);  

        %------------------------------------------------------------------
        %Debugging
        %figure; 
        %subplot(1,2,1); imagesc(small); axis image; colormap(gray(256));
        %subplot(1,2,2); imagesc(circular_edges); axis image; colormap(gray(256));
        %------------------------------------------------------------------
        
        if any(circular_edges(:))
            %if any edges are left, perform a hough circle search           
            rad_scores = hough_circle(circular_edges, [centres_x centres_y], min_rad, max_rad);
            
            %sort the hough totals to pick the most popular circles
            rad_scores_pos = sort(rad_scores(rad_scores > 0), 'descend');
            
            %for top three (and equal) most popular circles
            [pos_circ_cent pos_circ_rad] = find(rad_scores >= rad_scores_pos(3));
            
            %work out the corresponding centre and radii
            xc = zeros(length(pos_circ_cent), 1);
            yc = zeros(length(pos_circ_cent), 1);
            rad = zeros(length(pos_circ_cent), 1);
            for jj = 1:length(pos_circ_cent)
                xc(jj) = centres_x(pos_circ_cent(jj));
                yc(jj) = centres_y(pos_circ_cent(jj));
                rad(jj) = pos_circ_rad(jj) + min_rad - 1;
                %drawcircle(xc(jj),yc(jj),rad(jj),'g');
            end

            %Compute the final centre and radii as the average of these
            xc = mean(xc);
            yc = mean(yc);
            rad = mean(rad);
            
            %Draw the marker for the user
            manual_fig = figure(...
                'Name', 'Manually mark points',...
                'WindowStyle', 'normal',...
                'Position', [200 200 400 400],...
                'MenuBar', 'none');
            
            imagesc(small); axis image; colormap(gray(256));
            drawcircle(xc,yc,rad,'b');
            drawcross(xc,yc,10,'b');

            %Check it is ok        
            yes = 'yes';
            select_edge = 'No - manually select edge points';
            answer = questdlg(...
                'Is the selected circle OK?',...
                'Marker Point',yes, select_edge, yes);
            if strcmpi(answer, select_edge)
                reject = 1;
            end
            close(manual_fig);
        else
            reject = 1;
        end
    end        
    
    %if we have to manually select the points do that now
    while reject

        % display the contrast enhanced image
        manual_fig = figure(...
            'Name', 'Manually mark points',...
            'WindowStyle', 'normal',...
            'Position', [200 200 400 400],...
            'MenuBar', 'none');
        
        imagesc(small); axis image; colormap(gray(256));
        title('Click 3 or more points on the edge of the marker');
        xlabel('Press enter when finished');
        set(gca, 'xticklabel', [], 'yticklabel', []);
        
        %get the user to select points on the edge
        [x_select, y_select, P] = impixel; %#ok
        
        %Fit circle to these points
        [yc, xc, rad] = circfit(y_select,x_select);
        
        %Draw the circle on the image and check this is ok
        drawcircle(xc,yc,rad,'b');
        drawcross(xc,yc,10,'b');
        answer = questdlg(...
                'Is the selected circle OK?','Marker Point','Yes', 'No', 'Yes');
        if strcmpi(answer, 'yes')
            reject = 0;
        end
        close(manual_fig);
    end
    
    % translate back to original image coordinates
    xc = xc + x_start - 1;
    yc = yc + y_start - 1;

    % draw on main image the marker boundary and centre
    figure(marker_fig);
    set(marker_fig,'Name','Results of marker location');
    drawcircle(xc, yc, rad, 'r');
    drawcross(xc, yc, 10, 'r');

    