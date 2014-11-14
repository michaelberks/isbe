function [line_response, orientation, scale, normal_response] = gaussian_2nd_derivative_line(im, scales, g_width)
%GAUSSIAN_2ND_DERIVATIVE_LINE *Insert a one line summary here*
%   [line_response, orientation, scale] = gaussian_2nd_derivative_line(im, scales)
%
% Inputs:
%      im - *Insert description of input variable here*
%
%      scales - *Insert description of input variable here*
%
%
% Outputs:
%      line_response - *Insert description of input variable here*
%
%      orientation - *Insert description of input variable here*
%
%      scale - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 01-Dec-2010
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester

if nargin == 1
    %Get number of scales and image size from input responses
    responses = im; clear im;
    if iscell(responses)
        [r c] = size(responses{1}(:,:,1));
        num_scales = length(responses);
    else
        [r c num_scales] = size(responses(:,:,:,1));
    end
else
    %Set default filter width if none specified
    if nargin < 3
        g_width = 5;
    end
    
    %Get size of image then pad
    im = double(im);
    [r c] = size(im);
    padsize = round(g_width*max(scales));
    im = padarray(im, [padsize padsize], 'replicate');
    num_scales = length(scales);
end

%pre-allocate output arguments 
line_response = zeros(r, c);
orientation = zeros(r, c);
scale = zeros(r, c);
normal_response = zeros(r, c);

for i_scale = 1:num_scales
    
    if exist('responses', 'var')
        if iscell(responses)
            %Interpolate the responses back to the full pixel grid
            scaling = 2^(i_scale-1);
            offset = (i_scale-1) / (2^i_scale);
            rr = offset + (1:r)' / scaling;
            cc = offset + (1:c) / scaling;
            Ixy = interp2(responses{i_scale}(:,:,1), cc, rr, 'cubic');
            Ixx = interp2(responses{i_scale}(:,:,2), cc, rr, 'cubic');
            Iyy = interp2(responses{i_scale}(:,:,3), cc, rr, 'cubic');
        else
            Ixy = responses(:,:,i_scale,1);
            Ixx = responses(:,:,i_scale,2);
            Iyy = responses(:,:,i_scale,3);
        end
        sigma = i_scale;
    else
        %Make 2nd order directional filters
        sigma = scales(i_scale);
        [g,dg,ddg] = gaussian_filters_1d(sigma, round(g_width*sigma));

        %Filter the image
        % Filter the image with separable filters and remove padding
        Ixy = conv2(dg',dg,im,'same'); Ixy = -Ixy(padsize+(1:r), padsize+(1:c));
        Ixx = conv2(g',ddg,im,'same'); Ixx = Ixx(padsize+(1:r), padsize+(1:c));
        Iyy = conv2(ddg',g,im,'same'); Iyy = Iyy(padsize+(1:r), padsize+(1:c));
        
    end
    
    %theta_a = atan(2*Ixy ./ (Ixx-Iyy)) / 2; 
    theta_a = atan2(2*Ixy , (Ixx-Iyy)) / 2;
    theta_b = theta_a + pi/2;
    wo_theta_a = filter_output(Ixx, Iyy, Ixy, theta_a);
    wo_theta_b = filter_output(Ixx, Iyy, Ixy, theta_b);
    
    %if wo_theta_a at this scale so far is biggest, swap it into the main
    %line_response, orientation and scale matrices - note the theta are
    %normal to the line direction, so for theta_a we swap in theta_b and
    %vice-versa
    swap_idx = abs(wo_theta_a) > abs(line_response);
    line_response(swap_idx) = wo_theta_a(swap_idx);
    normal_response(swap_idx) = wo_theta_b(swap_idx);
    orientation(swap_idx) = theta_b(swap_idx);
    scale(swap_idx) = sigma;
    
    %Repeat for theta_b
    swap_idx = abs(wo_theta_b) > abs(line_response);
    line_response(swap_idx) = wo_theta_b(swap_idx);
    normal_response(swap_idx) = wo_theta_a(swap_idx);
    orientation(swap_idx) = theta_a(swap_idx);
    scale(swap_idx) = sigma;
    
end

%wrap orientations back into [-pi pi] range
orientation(orientation > pi) = orientation(orientation > pi) - 2*pi;

function wo_theta = filter_output(Ixx, Iyy, Ixy, theta)
%Compute the output of the filter for an arbitrary angle
cc = cos(theta).^2;
ss = sin(theta).^2;
s2 = sin(2*theta);
wo_theta = Ixx.*cc + Iyy.*ss + Ixy.*s2;

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Some old experimental code below, commented out but not ready for
%deletion...
% function [line_response, orientation, scale, minim, maxim] =...   
%   gaussian_2nd_derivative_line(im, scales)
% 
% f_debug = (nargin==0 && nargout==0);
% 
% if f_debug
% % 	datapath = 'A:\data\mammograms\2004_screening_processed\mass_roi';
% % 	load(fullfile(datapath,'045LML_roi.mat'));
% % 	im = bg(401:700,151:450);
% % 	scales = [1,2,4,8];
% 	
% 	datapath = 'U:\projects\nailfold\images\';
% 	im = mean(imread([datapath,'oneloop.png']),3);
% 	scales = [3];
% end
% 
% %pre-allocate output arguments
% [r c] = size(im);
% line_response = zeros(r, c);
% orientation = zeros(r, c);
% scale = zeros(r, c);
% 
% im = double(im);
% for i_scale = 1:length(scales)
%     
%     %Make 2nd order directional filters
%     [w0 w1 w2] = make_g2d_filt(scales(i_scale));
% 	
%     %Filter the image
% %     wo0 = imfilter(im, w0, 'same', 'replicate');
% %     wo1 = imfilter(im, w1, 'same', 'replicate');
% %     wo2 = imfilter(im, w2, 'same', 'replicate');
%     wo0 = imfilter(im, w0);
%     wo1 = imfilter(im, w1);
%     wo2 = imfilter(im, w2);
% 	
%     theta_a = (atan(sqrt(3) * (wo2 - wo1) ./ (wo1 + wo2 - 2*wo0))) / 2;
%     theta_b = theta_a + pi/2;
%     wo_theta_a = filter_output(wo0, wo1, wo2, theta_a);
%     wo_theta_b = filter_output(wo0, wo1, wo2, theta_b);
%     
% 	minim = min(wo_theta_a,wo_theta_b);
% 	maxim = max(wo_theta_a,wo_theta_b);
% 	
%     %if wo_theta_a at this scale so far is biggest, swap it into the main
%     %line_response, orientation and scale matrices - note the theta are
%     %normal to the line direction, so for theta_a we swap in theta_b and
%     %vice-versa
%     swap_idx = abs(wo_theta_a) > abs(line_response);
%     line_response(swap_idx) = wo_theta_a(swap_idx);
%     orientation(swap_idx) = theta_b(swap_idx);
%     scale(swap_idx) = scales(i_scale);
%     
%     %Repeat for theta_b
%     swap_idx = abs(wo_theta_b) > abs(line_response);
%     line_response(swap_idx) = wo_theta_b(swap_idx);
%     orientation(swap_idx) = theta_a(swap_idx);
%     scale(swap_idx) = scales(i_scale);
% end
% 
% %wrap orientations back into [-pi pi] range
% orientation(orientation > pi) = orientation(orientation > pi) - 2*pi;
% 
% if f_debug
% 	figure(1); clf; colormap(gray(256));
% 		imagesc(line_response);
% 		axis('image');
% 	clear; 
% end
% 
% function [w0 w1 w2] = make_g2d_filt(sigma)
% %aux function to make the gaussian 2nd derivative filters
% 
% width = round(5*sigma);
% ssq = sigma^2;
% 
% [x0 y0] = meshgrid(-width:width, -width:width);
% xy1 = [cos(pi/3) -sin(pi/3); sin(pi/3) cos(pi/3)] * [x0(:) y0(:)]';
% xy2 = [cos(2*pi/3) -sin(2*pi/3); sin(2*pi/3) cos(2*pi/3)] * [x0(:) y0(:)]';
% 
% x1 = reshape(xy1(1,:), 2*width+1, 2*width+1);
% y1 = reshape(xy1(2,:), 2*width+1, 2*width+1);
% x2 = reshape(xy2(1,:), 2*width+1, 2*width+1);
% y2 = reshape(xy2(2,:), 2*width+1, 2*width+1);
% 
% w0 = exp(-(x0.*x0)/(2*ssq)) .* exp(-(y0.*y0)/(2*ssq)) .* ((x0.*x0)/ssq - 1) * sqrt(2) / (sqrt(pi*ssq*ssq));
% w1 = exp(-(x1.*x1)/(2*ssq)) .* exp(-(y1.*y1)/(2*ssq)) .* ((x1.*x1)/ssq - 1) * sqrt(2) / (sqrt(pi*ssq*ssq));
% w2 = exp(-(x2.*x2)/(2*ssq)) .* exp(-(y2.*y2)/(2*ssq)) .* ((x2.*x2)/ssq - 1) * sqrt(2) / (sqrt(pi*ssq*ssq));
% 
% [sum(abs(w0(:))) sqrt(sum(w0(:).*w0(:))) max(w0(:))]
% % w0 = w0 / sum(w0(:));
% % w1 = w1 / sum(w1(:));
% % w2 = w2 / sum(w2(:));
% 
% function wo_theta = filter_output(wo0, wo1, wo2, theta)
% %Compute the output of the filter for an arbitrary angle
% c2 = cos(2*theta);
% s2 = sin(2*theta);
% wo_theta = ((1 + 2*c2) .* wo0 + ...
%             (1 - c2  + sqrt(3)*s2) .* wo1 + ...
%             (1 - c2  - sqrt(3)*s2) .* wo2) / 3;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------