function [xcbig,ycbig,R] = find_marker_legacy(IMAGE, xi, yi, reject, marker_fig)

% boxsize is side of square, centred on selected point, in which to search for marker 
    boxsize = 200;
    dimensions = size(IMAGE);
    
    option=0;
    
    while(option < 1.5)

        if(option==1)
            % select new box avoiding objects not part of marker
            disp('select top-left and bottom right corner of new search region, then press enter');
            new_box_fig = figure('Name', 'Select new search region', 'WindowStyle', 'normal');
            [newx, newy, P]=impixel(small); %#ok
            yi=round((newy(1)+newy(2))/2 + yi - boxsize/2 - 1);
            xi=round((newx(1)+newx(2))/2 + xi - boxsize/2 - 1);
            boxsize=round((newx(2)-newx(1))/2)*2;  % this ensures even boxsize
            close(new_box_fig);
       end

       % if box goes off edge of main image then nudge it in a bit
       % allow extra 10 pixels to avoid edges due to film edge
       if(xi-boxsize/2 < 1)
           xi = 1 + boxsize/2 + 10;
       elseif(xi+boxsize/2 > dimensions(2))
           xi = dimensions(2) - boxsize/2 - 10;
       end
       if(yi-boxsize/2 < 1)
           yi = 1 + boxsize/2 + 10;
       elseif(yi+boxsize/2 > dimensions(1))
           yi = dimensions(1) - boxsize/2 - 10;
       end
       
        % extract the region around this point (processing would be too slow to operate on whole image)
        small = IMAGE(yi-boxsize/2:yi+boxsize/2,xi-boxsize/2:xi+boxsize/2);
        
        % cast small image to double for manipulation
        dsmall = double(small);

        % contrast enhance the image segment
        minval = min(dsmall);
        minval = min(minval);
        maxval = max(dsmall);
        maxval = max(maxval);

        dsmall = dsmall-minval;
        dsmall = dsmall*255/(maxval-minval);

        % cast back to uint8 (0-255)
        small = uint8(dsmall);

        find_edge=1;
        thresholds = [0.1 0.5];
    
        while(find_edge == 1) 
            % apply canny edge detector to enhanced image
            cannysmall = edge(small,'canny',thresholds,2);
   
            %figure;
            %imshow(cannysmall);
        
            % if too few edge points are found then reduce thresholds and try again
            if(sum(sum(cannysmall))<40)
                find_edge=1;
                thresholds = thresholds/5;
            else
                find_edge=0;
            end
        end

        % make vector of points on the edge:
        [i,j] = find(cannysmall);
        
        % fit to circle
        [yc,xc,R] = circfit(i,j);

        % check quality of circle fit. 
        % Look for any edge points more than 10 pixels outside circle
        bad_points = sqrt(((j-xc).*(j-xc))+((i-yc).*(i-yc)))>R+10;
        if(sum(bad_points)>10)
            disp('number of bad points:');
            disp(sum(bad_points));
            select_edge=1;
        elseif(R>boxsize/2)
            disp('circle seems rather large!');
            disp(R);
            select_edge=1;
        elseif(xc+R>boxsize || xc-R<0 || yc+R>boxsize || yc-R<0)
            disp('circle goes outside box');
            disp([xc yc R]);
            select_edge=1;
        elseif(reject==1)
            disp('user requests manual override');
            select_edge=1;
        else
            select_edge=0;
            option=3;
        end
    
        if(select_edge)
            % draw the circle on the contrast enhanced image, and give option of accepting it.
            manual_fig = figure('Name', 'Check circle is OK',...
                'WindowStyle', 'normal', 'Position', [200 200 400 400]);
            imagesc(small); axis image; colormap(gray(256));
            drawcircle(xc,yc,R,'b');
            drawcross(xc,yc,10,'b');

            % either allow user to drag new box around circle which avoids any confusing objects,
            % or else to select some points on the circle edge.
            yes = 'yes';
            select_edge = 'No - select edge points';
            select_box = 'No - define new box';
            answer = questdlg(...
                'Is the selected circle OK?',...
                'Marker Point',yes, select_edge, select_box, yes);
            if strcmpi(answer, yes)
                option = 3;
            elseif strcmpi(answer, select_edge)
                option = 2;
            elseif strcmpi(answer, select_box)
                option = 1;
            end
%             disp('Select: "1" to define new box');
%             disp('        "2" to select edge points');
%             disp('        "3" to accept circle');
%             option = input(':');
            if option ~= 2
                close(manual_fig);
            end
        end
    
    end
    
    if(option==2)                   
        cont=0;
        while(cont==0)
            % allow user to select some points on the circle edge
            disp('select some points on the circle edge, then press enter');
            figure(manual_fig);
            set(manual_fig, 'Name','Select circle edge points');
            imagesc(small); axis image; colormap(gray(256));
            [x_select, y_select, P] = impixel; %#ok
            [yc, xc, R] = circfit(y_select,x_select);
            drawcircle(xc,yc,R,'b');
            drawcross(xc,yc,10,'b');
            answer = questdlg(...
                'Is the selected circle OK?','Marker Point','Yes', 'No', 'Yes');
            if strcmpi(answer, yes)
                cont = 1;
            end
        end
        close(manual_fig);
    end
    
    % translate back to original image coordinates
    xcbig = xc + xi - boxsize/2 - 1;
    ycbig = yc + yi - boxsize/2 - 1 ;
    
    
    %disp('Centre of circle in small image is at: ')
    %disp([xc,yc])
    %disp('Centre of circle in main image is at: ')
    %disp([xcbig,ycbig])
    %pixval

    % draw on main image the marker boundary and centre
    figure(marker_fig);
    set(marker_fig,'Name','Results of marker location');
    drawcircle(xcbig,ycbig,R,'r');
    drawcross(xcbig,ycbig,10,'r');

    