function [out_var, scaling] = load_roi8(filepath, row, col, win_size, varargin)
%LOAD_UINT8 load a variable that has been save using SAVE_UINT8 - this will
%undo the scaling and convert the variable back into double format
%   [out_var, scaling] = load_uint8(filepath)
%
% Inputs:
%      filepath - of variable saved using SAVE_UINT8
%
%
% Outputs:
%      out_var - output variable scaled and converted into double format
%
%      scaling - scaling stored with saved variable
%
%
% Example:
%
% Notes:
%
% See also: SAVE_UINT8
%
% Created: 27-Oct-2010
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester

% check for values that are identical to that mapped from zero
check_zero = false;
if nargin>4, check_zero = varargin{1}; end

s = load(filepath);

if ~all(isfield(s, {'out_var', 'scaling'}))
    warning('load_uint8:missing_variables', 'Variable out_var and scaling not present in the saved file. Are you sure it was save using SAVE_UINT8');
    names = fieldnames(s);
    if length(names) == 1
        out_var = getfield(s, names{1}); %#ok
    else
        out_var = s;
    end
    return;
end

scaling = s.scaling;

do_complex = isfield(scaling, 'complex');

if do_complex
    if isfield(s, 'zero_mask')
        dims = size(s.zero_mask);
        real_mask = sample_window(s.zero_mask, win_size, row, col, 0);
        imag_mask = sample_window(s.zero_mask, win_size, row + dims(1)/2, col, 0);
        s = rmfield(s, 'zero_mask');
        
        real_out = zeros(win_size);
        imag_out = zeros(win_size);
        
        real_out(real_mask) = 
    
else
    dims = size(out_var);
    real_part = out_var(1:dims(1)/2,:);
    imag_part = out_var(1+dims(1)/2:dims(1),:);
    
    out_var = reshape(complex(real_part, imag_part), [dims(1)/2 dims(2:end)]);
end

if isfield(s, 'zero_mask');
    out_var = double(s.zero_mask);
    out_var(s.zero_mask) = (scaling.m * double(s.out_var) / 255) + scaling.c;
else
    out_var = (scaling.m * double(s.out_var) / 255) + scaling.c;
end

% find values that look they might have been zero before compression and
% set to zero rather than some small number that has been corrupted by
% rounding error
if (check_zero)
	zero_val = uint8(255 * (0-scaling.c)/scaling.m);
	out_var(s.out_var == zero_val) = 0;
end


