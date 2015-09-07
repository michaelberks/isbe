function [c_array] = read_complex_txt(filename)
%READ_COMPLEX_TXT read in a 2D complex array from ascii text format
%   [] = read_complex_txt(complex_mat, fname)
%
% Inputs:
%      filename - of text to be read
%
%
% Outputs: c_array - 2d complex array read from the file
%
% Example:
%
% Notes: Each element C=a+ib must be written in the form (a, b). The elements are space delimited, each row starting a new line
% For some reason, text written by our C++ implementaton (the main reason
% for the existence of this function) doesn't pick up the end-of-line
% characters so we process these manually. This does mean the function may
% not work with text files written from other formats
% See also:
%
% Created: 23-Apr-2015
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
fid = fopen(filename,'r');
frewind(fid);
s = textscan(fid, '%s', 'endofline', '\n', 'delimiter', {' ', ', ', '\n'});
fclose(fid);
s = s{1};

%Merge real and imaginary parts of the non-zero elements together
non_zeros_r = strncmpi(s, '(', 1);
if s{find(non_zeros_r,1)}(end) ~= ')'
    non_zeros_i = [false; non_zeros_r(1:end-1,:)];
    s(non_zeros_i) = strcat(s(non_zeros_r), ',', s(non_zeros_i));
    s(non_zeros_r) = [];
end

%Find the end-of-lines, for some reason the initial textscan won't
%recognise these. Discard and check the expected row and column sizes match
eol = ~strncmpi(s, '(', 1) & ~strncmpi(s, '0', 1);
if any(eol)
    col_size = find(eol,1)-1;
    s(eol) = [];
    row_size = length(s) / col_size;
    %if row_size ~= sum(eol)
    %    error('End of line paramaters incoorectly read');
    %end
    %Otherwise resize s (note s will be a column vector, in row order, so
    %need to transpose)
    s = reshape(s, col_size, row_size)';
end

%Can now read in each complex element
c_array = zeros(size(s));
for i_r = 1:size(s,1)
    for i_c = 1:size(s,2)
        if ~strncmpi(s{i_r,i_c}, '0', 1)           
            c_num = textscan(s{i_r,i_c}, '(%n,%n)');
            c_array(i_r,i_c) = complex(c_num{1}, c_num{2});
        end
    end
end


