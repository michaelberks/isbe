function [ markers ] = single_separation_check(markers)
%MARKER_CHECK : Given input markers, function performs discrimination against these markers and outputs the results
%Assumes we have at most 1 marker in each region
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
markers_a = markers(:,3)==1;
markers_b = markers(:,3)==2;

num_a = sum(markers_a); 
num_b = sum(markers_a);

%Compare separations
if num_a && num_b %if there are markers on either side of the image
    
    regions_a = unique(markers(markers_a,4));
    regions_b = unique(markers(markers_b,4));
    
    if length(regions_a)<2 || length(regions_b)<2
        return;
    end
    
    for i_r = regions_a(1:end-1)'
        
        %Need i_r and i_r + 1 in both A and B
        if ismember(i_r+1, regions_a) && ismember(i_r, regions_b) && ismember(i_r+1, regions_b)
            %By definition, each of the below index vectors should have
            %exactly one true value
            idx_a = (markers(:,4) == i_r) & markers_a; 
            idx_a1 = (markers(:,4) == i_r+1) & markers_a;
            idx_b = (markers(:,4) == i_r) & markers_b;
            idx_b1 = (markers(:,4) == i_r+1) & markers_b;
            
            idx_all = idx_a | idx_a1 | idx_b | idx_b1;
            
            %Increment tests performed
            markers(idx_all,6) =  markers(idx_all,6) + 1;
                
            %Calculate separations
            a_sep = eucl_distance(...
                markers(idx_a,1), markers(idx_a,2),...
                markers(idx_a1,1), markers(idx_a1,2)); 
            b_sep = eucl_distance(...
                markers(idx_b,1), markers(idx_b,2),...
                markers(idx_b1,1), markers(idx_b1,2)); 
            
            %Compute ratio and check against bounds
            sep_ratio = a_sep/b_sep;
            if sep_ratio >=lowerbound && sep_ratio<=upperbound 
                markers(idx_all,5) = markers(idx_all,5) + 1;
            end
        end
    end
end