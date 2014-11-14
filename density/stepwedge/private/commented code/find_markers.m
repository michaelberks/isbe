function [ image,markers,x_s,y_s,delete_range,region ] = find_markers( oimage , plot_results ,extra_plots)
%FIND_MARKERS : Highest order function used to find markers given the original image
% 
%
% Inputs:
%			oimage				image that you want to find markers in [2D array]
%			plot_results 		1 if you want to display found markers on image, 0 if not [#]
%			extra_plots			1 if you want to display extra intermediate images used for analysis of result, 0 if not [#]
%
% Outputs:
%			image				the original image that has now been cut down in size by removing the white frame around its edge [2D array]
%			markers 			x-y position of the markers found by the software [2D array] 
%								x-position/y-position/strip number/region ID number output
%			x_s 				gives the x-split of the edge strips taken for analysis within a single image (varies with image size) [1D array]
%			y_s 				gives the y-split of the edge strips taken for analysis within a single image (varies with image size) [1D array]
%			delete_range		number of pixel colomns or rows removed from the image edges during "remove_white_frame.m"
%			region 				which set of small regions contains the most found markers [#]
%
% Example:
%
% Notes:	I think it would be best to make the white frame around the image a mask rather than removing it. As the code stands at the moment
%			the marker co-ordinates are in the frame of the cut image and not the original. 
%			
% See also:
%
% Created: 18-06-2013
% Author: Euan Allen
% Email : euan.allen@gmail.com



edge_sample = 600; %size of image edge strip that the code searches for the markers in
zero_region = 50; %'zero region' generated around found peaks

if length(oimage)>6000; %repeat for small images, not large ones (due to regions)
    R=1; 
else
    R=2;
end

for i = 1:R %two region definitions
    %Prepare image and find the allow_map--------------------------------------
    [r_map] = region_definitions(oimage,i); %Define regions
    [image, r_map,delete_range] = remove_white_frame(oimage,r_map); %Remove white frame around edge
    [image, a_map] = anon_white(image); %Make anon white and supply anon map
    allow_map = r_map.*a_map;%Allow map eliminates if some of a region goes in anon area
    [Y,X] = size(image);
    clear r_map
    %--------------------------------------------------------------------------
    
    %Region splits-------------------------------------------------------------
    if length(image)>6000 %large image
        x_s = [1, X, 1, X]; %take 0->X for both splits
        y_s = [1,edge_sample , Y-edge_sample+1, Y]; %edge co-ords of split
        number_of_peaks = 4; %initial peaks to consider within an image strip
    else%small image
        x_s = [1,edge_sample,X-edge_sample+1,X];
        y_s = [1,Y 1, Y];
        number_of_peaks = 3;
    end
     %--------------------------------------------------------------------------
	 
    %Find markers in split regions----------------------------------------------
    [markers_raw] = find_algorithm(image,x_s,y_s,zero_region,number_of_peaks,allow_map,a_map,extra_plots);
	 %--------------------------------------------------------------------------
	
	%Discriminate if possible---------------------------------------------------
    if markers_raw(1,:) ~= 0 %if markers have been found
        temp_markers = marker_check(markers_raw);%separation check (discrimination)
    else
        display('No markers found')
        temp_markers = markers_raw;
    end
	%--------------------------------------------------------------------------
	
    marker(i).markers = temp_markers; %required if repeating for both region definitions for small images
    marker(i).numbermarkers = size(temp_markers,1);
    clear temp_markers
end

%Find which regional set is the likely correct one (based on number of markers found in each)
if R==1 %large image (only one region definition)
    region = 1;
    markers = marker(1).markers;
elseif marker(1).numbermarkers>marker(R).numbermarkers %if region 1 has more markers than 2
    region = 1; %which region is the likely current one
    markers = marker(1).markers;	%set output as markers found in that region
elseif marker(1).numbermarkers==marker(R).numbermarkers%if region 1 has the same number of markers than 2
    region = 1;
    display('Both regions have same number of found markers')
    markers = [marker(1).markers; marker(R).markers]; 
else	%if region 2 has more markers than 1
    region = 2;
    markers = marker(R).markers;
end
%--------------------------------------------------------------------------

%Plot results if require by function input
if plot_results
    close
    figure;imshow(image)
    hold on
    for C = 1:size(markers,1)
        if markers(C,5) == 0 %failed all separations
            plot(markers(C,1),markers(C,2),'r.','MarkerSize',10) %plot as red dot
        elseif markers(C,5) ~= markers(C,6) %failed some separations but not all
            plot(markers(C,1),markers(C,2),'b.','MarkerSize',10); %plot as blue dot
        else %passed all discriminations
            plot(markers(C,1),markers(C,2),'c.','MarkerSize',10) %plot as cyan dot
        end
    end
    pause
end
%--------------------------------------------------------------------------

end

