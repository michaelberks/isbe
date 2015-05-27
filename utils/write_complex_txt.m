function [] = write_complex_txt(complex_mat, fname)
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
fid1 = fopen(fname, 'wt');
for i_r = 1:size(complex_mat, 1)
    for i_c = 1:size(complex_mat,2)
        if abs(complex_mat(i_r,i_c)) > 1e-4
            fprintf(fid1,'(%.4f, %.4f) ', real(complex_mat(i_r,i_c)), imag(complex_mat(i_r,i_c)));
        else
            fprintf(fid1, '0 ');
        end
    end
    fprintf(fid1,'\n');
end
fclose(fid1);


