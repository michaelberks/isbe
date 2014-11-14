function [markers, mValues] = template_matching(image,template,zero_range,top_sample_number,mask,display)
%TEMPLATE_MATCHING : Performs template matching algorithm to an image
% 
%
% Inputs:
%			image				image to perform template match on
%			template 			template to use during template matching
%			zero_range			region to zero around a found peak
%			top_sample_number	number of peaks wanted to output
%			mask	 			mask to show regions where a peak may exist (not anon region)
%			display				set to 1 to display results, 0 otherwise
%
% Outputs:
%			markers 			x-y position of the markers found by the software [2D array] 
%			mValues 			corresponding cross correlation value for the peaks found [# (between -1 and 1)]
%
% Example: [markers, mValues] = template_matching(image,'patch3.mat',200,6,1)
%
% Notes:	Function requires access to 'patch3.mat'.
%			
% See also:
%
% Created: 18-06-2013
% Author: Euan Allen
% Email : euan.allen@gmail.com


%Sets the image into correct format (if in RBG)
if length(size(image)) >= 3                %If not already greyscaled
    Io = image(:,:,1);                  %Convert to greyscale image
elseif length(size(image)) == 2             %If greyscaled
    Io = image;
else                                        %If not in right format
    error('Image not correct format once read in')
end
%-------------------------------------------

%Add `white' frame (`white' --> 2^12 to shift the correlation
%frame to outside of the original image
f=200; %frame size
[Ao Bo] = size(Io);
I = 65536*ones(Ao+2*f,Bo+2*f); %create a base of x*ones with a 200 frame around edge.
I(f+1:f+Ao,f+1:f+Bo) = Io(:,:);
%-------------------------------------------

%Cross correlate
template_patch = load(template,'patch'); %load template
corr_map = normxcorr2(template_patch.patch,I); %cross correlate template with image
%-------------------------------------------

%Remove corr_map padding (MATLAB added zeros)
[A B] = size(corr_map);
pad = 0.5*([A B]-size(I));
corr_map_NP = corr_map(pad(1)+1:A-pad(1),pad(2)+1:B-pad(2));
%-------------------------------------------

%Remove self added frame
corr_map_NP = corr_map_NP(f+1:f+Ao,f+1:f+Bo);
%-------------------------------------------

%Find Maxima
[tempmarkers tempmValues] = local_image_maxima(corr_map_NP,zero_range,mask);
%-------------------------------------------

if length(tempmarkers)>top_sample_number %if found more than wanted
    markers = tempmarkers(1:top_sample_number,:); %take the amount you do want starting from the largest peak 
    mValues = tempmValues(1:top_sample_number);
else %take all found
    markers = tempmarkers;
    mValues = tempmValues;
end

if display %if user wants to display results
    figure;imshow(Io)
    hold on
    for C = 1:size(markers,1)
        plot(markers(C,1),markers(C,2),'r.','MarkerSize',10)
        text(markers(C,1)+20,markers(C,2)+20, num2str(C),'Color',[1 1 1],'BackgroundColor',[0 0 0])
    end    
end

end


