function [] = write_array_txt(in_array, fname)
%WRITE_COMPLEX_TXT write out a 2D array in ascii text format
%   [] = write_array_txt(complex_mat, fname)
%
% Inputs:
%      in_array - 2d array to be written
%
%      fname - filename to be written to
%
%
% Outputs: none
%
% Example:
%
% Notes: Each element 4 digits of
% precision. If, an element is less than 1e-4 it is written as 0 to save
% space
%
% See also:
%
% Created: 23-Apr-2015
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
fid1 = fopen(fname, 'wt');
for i_r = 1:size(in_array, 1)
    for i_c = 1:size(in_array,2)
        if abs(in_array(i_r,i_c)) > 1e-4
            fprintf(fid1,'%.4f ', in_array(i_r,i_c));
        else
            fprintf(fid1, '0 ');
        end
    end
    fprintf(fid1,'\n');
end
fclose(fid1);


