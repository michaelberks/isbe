function [BW] = icp_hysterisis(icp_band, ilp_band, mag_thresh, phase_thresh, plot_flag)
%HYSTERISIS *Insert a one line summary here*
%   [BW] = hysterisis(magnitudes,thresh)
%
% Inputs:
%      magnitudes - *Insert description of input variable here*
%
%      mag_thresh - *Insert description of input variable here*
%
%      phase_thresh -
%
%
% Outputs:
%      BW- *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 07-Nov-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester 
if nargin < 5
    plot_flag = 0;
end
if nargin < 4
    phase_thresh = [];
end
if isempty(phase_thresh)
    phase_thresh = [-inf -inf];
end

magnitudes = abs(icp_band);
icp_phase = angle(icp_band);
ilp_phase = angle(ilp_band);
%magnitudes = magnitudes / max(magnitudes(:));

[m n] = size(magnitudes);

%Do NMS - get normal orientations
ori_c = sin(icp_phase);
ori_s = cos(icp_phase);

%compute (x, y) coordinates along +/- unit normal vectors
xx = repmat(1:n, m, 1);
yy = repmat((1:m)', 1, n);

x_interp1 = xx + ori_c;
x_interp2 = xx - ori_c;

y_interp1 = yy + ori_s;
y_interp2 = yy - ori_s;

%linearly interpolate values at normal coordinates
z1 = interp2(magnitudes, x_interp1, y_interp1);
z2 = interp2(magnitudes, x_interp2, y_interp2);

discard = (magnitudes <= z1) | (magnitudes <= z2);
nms_mag = magnitudes;
nms_mag(discard) = 0;
    
%Now do the edge hysterisis
%Compute the strong edges and the weak edges
weak_lines = (nms_mag > mag_thresh(1)) & (ilp_phase > phase_thresh(1));
strong_lines = (nms_mag > mag_thresh(2)) & (ilp_phase > phase_thresh(2));

%Now find weak edges that are 8-connected to a strong edge
if any(strong_lines(:))
    [rstrong cstrong] = find(strong_lines);
    combined_lines = bwselect(weak_lines, cstrong, rstrong, 8);
else
    combined_lines = weak_lines;
end

if plot_flag
    figure; 
    subplot(2,2,1); imagesc(nms_mag); axis image; colormap(gray(256));
    subplot(2,2,2); imagesc(weak_lines); axis image; colormap(gray(256));
    subplot(2,2,3); imagesc(strong_lines); axis image; colormap(gray(256));
    subplot(2,2,4); imagesc(combined_lines); axis image; colormap(gray(256));

%     icp_mag = histeq(magnitudes / max(magnitudes(:)), 1 ./ (1:256));
%     figure; 
%     image(complex2rgb(icp_mag.*exp(i*2*icp_phase))); axis image; hold on;
%     for tx = 1:n; 
%         for ty = 1:m;
%             if combined_lines(ty,tx)
%                 plot([tx x_interp2(ty,tx)], [ty y_interp2(ty,tx)], 'g');
%                 plot([tx x_interp1(ty,tx)], [ty y_interp1(ty,tx)], 'y');
%             end
%         end;
%     end;
end

%finally, thin the combined edges
%BW = bwmorph(combined_lines, 'skel', inf);
%BW = combined_lines;
BW = bwmorph(combined_lines, 'thin', 1);