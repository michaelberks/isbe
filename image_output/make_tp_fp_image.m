function [tp_fp_image] = make_tp_fp_image(true_label, false_label, sample_label)
%MAKE_TP_FP_IMAGE *Insert a one line summary here*
%   [tp_fp_image] = make_tp_fp_image(true_label, sample_label)
%
% Inputs:
%      true_label - *Insert description of input variable here*
%
%      sample_label - *Insert description of input variable here*
%
%
% Outputs:
%      tp_fp_image - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 21-Jan-2011
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester

%Pre-allocate r g b channels (will be combined at end of function)
[r c] = size(true_label);
tp_fp_r = ones(r,c);
tp_fp_g = ones(r,c);
tp_fp_b = ones(r,c);

%1) Set TN (i.e. correctly identified background) whilte - don't need to
%change anything

%2) Set TP (i.e correctly identified foreground) red
tp = true_label & sample_label;
tp_fp_r(tp) = 1;
tp_fp_g(tp) = 0;
tp_fp_b(tp) = 0;

%3) Set FP (background mis-classified as foreground) green
fp = false_label & sample_label;
tp_fp_r(fp) = 0;
tp_fp_g(fp) = 1;
tp_fp_b(fp) = 0;

%4) Set FN (foreground mis-classified as background) blue
fn = true_label & ~sample_label;
tp_fp_r(fn) = 0;
tp_fp_g(fn) = 0;
tp_fp_b(fn) = 1;


%Combine channels into single output
tp_fp_image = cat(3, tp_fp_r, tp_fp_g, tp_fp_b);

