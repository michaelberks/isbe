function [] = write_vessel_contour_txt(vessel_contour, fname, scale_factor, reflect)
%WRITE_COMPLEX_TXT write out a 2D complex array in ascii text format
%   [] = write_complex_txt(complex_mat, fname)
%
% Inputs:
%      complex_mat - 2d complex array to be written
%
%      fname - filename to be written to
%
%
% Outputs: none
%
% Example:
%
% Notes: Each element C=a+ib will be written in the form (a, b) with 8 digits of
% precision. The elements are space delimited, each row starting a new line
%
% See also:
%
% Created: 23-Apr-2015
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester

if ~exist('scale_factor', 'var')
    scale_factor = 1;
end
if ~exist('reflect', 'var')
    reflect = false;
end

fid1 = fopen(fname, 'wt');

num_pts = size(vessel_contour.inner_edge,1);
rows = ceil(vessel_contour.rows*scale_factor);
cols = ceil(vessel_contour.cols*scale_factor);

vessel_contour.vessel_centre = vessel_contour.vessel_centre*scale_factor;
vessel_contour.inner_edge = vessel_contour.inner_edge*scale_factor;
vessel_contour.outer_edge = vessel_contour.outer_edge*scale_factor;

if reflect
    vessel_contour.vessel_centre(:,1) = cols - vessel_contour.vessel_centre(:,1) + 1;
    vessel_contour.inner_edge(:,1) = cols - vessel_contour.inner_edge(:,1) + 1;
    vessel_contour.outer_edge(:,1) = cols - vessel_contour.outer_edge(:,1) + 1;
end

apex_idx = vessel_contour.apex_idx;
apex = vessel_contour.vessel_centre(apex_idx,:);

fprintf(fid1, 'apex: %.4f %.4f\n', apex(:,1), apex(:,2));
fprintf(fid1, 'apex_idx: %d\n', apex_idx);
fprintf(fid1, 'rows: %d\n', rows);
fprintf(fid1, 'cols: %d\n', cols);

fprintf(fid1, 'vessel_centre: {\n');
for i_pt = 1:num_pts
    
    fprintf(fid1,'%.4f %.4f ', vessel_contour.vessel_centre(i_pt,1), vessel_contour.vessel_centre(i_pt,2));
    fprintf(fid1,'\n');
end
fprintf(fid1, '} //vessel_centre\n');

fprintf(fid1, 'inner_edge: {\n');
for i_pt = 1:num_pts
    
    fprintf(fid1,'%.4f %.4f ', vessel_contour.inner_edge(i_pt,1), vessel_contour.inner_edge(i_pt,2));
    fprintf(fid1,'\n');
end
fprintf(fid1, '} //inner_edge\n');

fprintf(fid1, 'outer_edge: {\n');
for i_pt = 1:num_pts
    
    fprintf(fid1,'%.4f %.4f ', vessel_contour.outer_edge(i_pt,1), vessel_contour.outer_edge(i_pt,2));
    fprintf(fid1,'\n');
end
fprintf(fid1, '} //outer_edge\n');

fclose(fid1);


