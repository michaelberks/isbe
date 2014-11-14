function [ markers_checked ] = marker_check(markers)
%MARKER_CHECK : Given input markers, function performs discrimination against these markers and outputs the results
% 
%
% Inputs:
%			markers 			x-y position of the markers found by the software [2D array] 
%								x-position/y-position/strip number/region ID number output XXXXXXXXXXXXXX
%
% Outputs:
%			markers_checked		x-y position of the markers found by the software [2D array] 
%								x-position/y-position/strip number/region ID number output/# of discrim. tests passed/# of discrim. tests performed
%
% Example:
%
% Notes:	
%			
% See also:
%
% Created: 18-06-2013
% Author: Euan Allen
% Email : euan.allen@gmail.com


lowerbound = 0.94; 
upperbound = 1.06; %upper and lower bands of the separation ratios (obtained via previous mammogram data)
markers(:,5:6) = 0; %add extra columns to markers

%Split a and b -type markers into two variables (so you can compare the separations)
markers_a = markers(markers(:,3)==1,:);
markers_b = markers(markers(:,3)==2,:);

num_a = size(markers_a,1); 
num_b = size(markers_a,1); 



%Compare separations
if num_a && num_b %if there are markers on either side of the image
    
    regions_a = unique(markers(:,4));
    regions_b = unique(markers(:,4));
    
    for i_ra = regions_a(1:end-1)'
        
        %Need i_ra and i_ra + 1 in both A and B
        if ismember(i_ra+1, regions_a) && ismember(i_ra, regions_b) && ismember(i_ra+1, regions_b)
            idx_a = markers_a(:,4) == i_ra;
            idx_a1 = markers_a(:,4) == i_ra+1;
            idx_b = markers_b(:,4) == i_ra;
            idx_b1 = markers_b(:,4) == i_ra+1;
        
    %sort by regional number
    [markers_a_sorted] = sortrows(markers_a, 4); 
    [markers_b_sorted] = sortrows(markers_b, 4);
    
    
    for i_a = 1:num_a - 1
        for i_b = 1:num_b-1
            
            %opposite pair with separations
            if markers_a_sorted(i_a,4) == markers_b_sorted(i_b,4)...
                    && markers_a_sorted(i_a+1,4)==markers_b_sorted(i_b+1,4) 
                
				%Calculate separations
                a_sep = eucl_distance(...
                    markers_a_sorted(i_a,1),markers_a_sorted(i_a,2),...
                    markers_a_sorted(i_a+1,1),markers_a_sorted(i_a+1,2)); 
                b_sep = eucl_distance(...
                    markers_b_sorted(i_b,1),markers_b_sorted(i_b,2),...
                    markers_b_sorted(i_b+1,1),markers_b_sorted(i_b+1,2));
                
                %Add 1 to total no. of tests
                markers_a_sorted(i_a,6)= markers_a_sorted(i_a,6)+1;
                markers_a_sorted(i_a+1,6)= markers_a_sorted(i_a+1,6)+1; 
                markers_b_sorted(i_b,6)= markers_b_sorted(i_b,6)+1;
                markers_b_sorted(i_b+1,6)= markers_b_sorted(i_b+1,6)+1;
                
                %If ratio of the two is between the condition (i.e. a pass)
                
                %add 1 to total no. of tests passed
                if a_sep/b_sep >=lowerbound && a_sep/b_sep<=upperbound 
                    markers_a_sorted(i_a,5)= markers_a_sorted(i_a,5)+1;
                    markers_a_sorted(i_a+1,5)= markers_a_sorted(i_a+1,5)+1; 
                    markers_b_sorted(i_b,5)= markers_b_sorted(i_b,5)+1;
                    markers_b_sorted(i_b+1,5)= markers_b_sorted(i_b+1,5)+1;
                end
            end
        end
    end
    markers_checked = [markers_a_sorted; markers_b_sorted]; %collect markers together

elseif num_a && ~num_b
    markers_checked = markers_a;

elseif ~num_a && num_b
    markers_checked = markers_b;
end


