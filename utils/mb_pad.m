function [array_out] = mb_pad(array_in, pad_size0, pad_val, where_str)
%MB_PAD *Insert a one line summary here*
%   [array_out] = mb_pad(in_array, pad_size, pad_val)
%
% Inputs:
%      in_array - *Insert description of input variable here*
%
%      pad_size - *Insert description of input variable here*
%
%      pad_val - *Insert description of input variable here*
%
%
% Outputs:
%      array_out - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 25-Jul-2011
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

if (nargin==0 && nargout==0), test_script(); return; end

% Default parameter values
if ~exist('pad_val','var') || isempty(pad_val), pad_val = 0; end
if ~exist('where_str','var') || isempty(where_str), where_str = 'both'; end

if ~strcmp(where_str, 'pre') && ...
   ~strcmp(where_str, 'post') && ...
   ~strcmp(where_str, 'both')
    error('Parameter 5 not valid: must be ''pre'', ''post'' or ''both''');
end

% Determine size of padded array
sz_in = size(array_in);
nDims = length(sz_in);

% Treat scalar pad_size as a padding on every dimension.
% We probably shouldn't do this - a user might actually want to pad on
% the first dimension only. Consider removing this at some point.
if (length(pad_size0) == 1)
    pad_size0 = pad_size0 * ones(1,nDims);
end

if (length(pad_size0) >= nDims)
    pad_size = pad_size0(1:nDims);
else
    pad_size = [pad_size0 zeros(1, nDims-length(pad_size0))];
end

% Ensure pad_size0 is an integer.
pad_size0 = ceil(pad_size0);

% Build up a mask via kronecker products, starting with scalar 'true'.
mask = true;
mask_sz = ones(1,nDims);

% Define initial size of output (may be modified).
sz_out = sz_in;

% If any elements of pad_size are negative then trim these dimensions
for i = 1:nDims
    % Number by which to pad/trim
    n = pad_size(i);
    
    f_pad = (0 < n);
    f_trim = (n < 0);
    
    if (f_trim)
        pad_size(i) = 0; % Ignore padding from now on.
    end

    % Determine the output size
    % Note that the matrix is reordered at every iteration so we always 
    % work on the first dimension.
    
    switch (where_str)
        case {'pre','post'},    sz_out(1) = max(0, sz_in(1)+n);
        case {'both'},          sz_out(1) = max(0, sz_in(1)+2*n);
    end
    
    % Set the mask for this dimension
    if f_pad
        mask_d = false(1, sz_out(1));

        switch (where_str)
            case {'pre'},   mask_d(1+n:end  ) = true;
            case {'post'},  mask_d(1  :end-n) = true;
            case {'both'},  mask_d(1+n:end-n) = true;
        end
    elseif f_trim
        mask_d = true(1, sz_out(1));

        % Trim values from the input matrix
        n = abs(n);
        sz_in(1) = sz_out(1);
        switch (where_str)
            case {'pre'},   array_in = reshape(array_in(1+n:end  , :), sz_in);
            case {'post'},  array_in = reshape(array_in(1  :end-n, :), sz_in);
            case {'both'},  array_in = reshape(array_in(1+n:end-n, :), sz_in);
        end
    else
        % Leave unchanged
        mask_d = true(1, sz_out(1));
    end

    % Build up mask one dimension at a time.
    mask_sz(i) = length(mask_d);
    mask = reshape(kron(mask_d(:)', mask), mask_sz);
    
    % Reorder to make next dimension first.
    array_in = shiftdim(array_in, 1);
    sz_in = circshift(sz_in, [0, -1]);
    sz_out = circshift(sz_out, [0, -1]);
end

mask = logical(mask);

% Return if no more padding to apply.
if all(pad_size == 0)
    array_out = array_in;
    return;
end

% Copy values across
array_out = zeros(sz_out);
array_out(mask) = array_in;

% Numeric pad value
if isnumeric(pad_val)
    if (pad_val ~= 0)
        array_out(~mask) = pad_val;
    end
    
    return
end

% Set repeat flag for symmetric padding.
sym_repeat = false;
if strcmp(pad_val, 'symmetric_repeat') || ...
   strcmp(pad_val, 'symmetric-repeat') % backward compatibility
    pad_val = 'symmetric';
    sym_repeat = true;
end

% Semantic pad definition (see padarray): 
%   'replicate'     Repeats border elements of A
%   'circular'      Pads with circular repetition of elements
%   'symmetric'     Pads array with mirror reflections of itself.
switch pad_val
    case 'replicate',
        for i = 1:nDims
            n = pad_size(i);
            if (strcmp(where_str,'pre') || strcmp(where_str,'both'))
                % Replicate beginning
                ind0 = 1 + n;
                inds = ind0 + zeros(1,n);
                array_out(1:n, :) = array_out(inds, :);
            end
            if (strcmp(where_str,'post') || strcmp(where_str,'both'))
                % Replicate end
                ind0 = sz_out(i) - n;
                inds = ind0 + zeros(1,n);
                array_out(end-n+1:end, :) = array_out(inds, :);
            end
            
            array_out = shiftdim(array_out, 1);
        end
        
    case 'circular',
        for i = 1:nDims
            n = pad_size(i);
            if (strcmp(where_str,'pre') || strcmp(where_str,'both'))
                % Replicate beginning
                ind0 = 1 + n;
                inds = ind0 + mod(-n:-1, sz_in(i));
                array_out(1:n, :) = array_out(inds, :);
            end
            if (strcmp(where_str,'post') || strcmp(where_str,'both'))
                % Replicate end
                ind0 = sz_out(i) - n - sz_in(i) + 1;
                inds = ind0 + mod(0:n-1, sz_in(i));
                array_out(end-n+1:end, :) = array_out(inds, :);
            end
            
            array_out = shiftdim(array_out, 1);
        end

    case 'symmetric',
        for i = 1:nDims
            n = pad_size(i);
            if (strcmp(where_str,'pre') || strcmp(where_str,'both'))
                % Replicate beginning
                ind0 = n;
                inds = reflect_inds(n, sz_in(i), sym_repeat);
                inds = ind0 + inds(end:-1:1); % reverse
                array_out(1:n, :) = array_out(inds, :);
            end
            if (strcmp(where_str,'post') || strcmp(where_str,'both'))
                % Replicate end
                ind0 = sz_out(i) - n + 1;
                inds = ind0 - reflect_inds(n, sz_in(i), sym_repeat);
                array_out(end-n+1:end, :) = array_out(inds, :);
            end
            
            array_out = shiftdim(array_out, 1);
        end

    otherwise,
        error(['Pad style ', pad_val, ' not recognized']);
end


function inds = reflect_inds(n, k, repeat)
% Generate sequence of n indices:
% repeat = true:  inds = 1, 2, ..., k-1, k, k, k-1, ..., 2, 1, 1, 2, ...
% repeat = false: inds = 2, 3, ..., k-1, k, k-1, ..., 2, 1, 2, ...

inds = [];

if (n == 0), return; end

if (repeat)
    inds = [1:k, k:-1:1];
else
    inds = [2:k, k-1:-1:1]; % Trim the first value
end

if isempty(inds), return; end

n_reps = ceil(n / length(inds));
inds = repmat(inds, [1,n_reps]);
inds = inds(1:n);


function test_script()
clc;

where_str = 'post';

m = reshape(1:25, [5,5]);
m2 = mb_pad(m, 0, 'symmetric', where_str)
m2 = mb_pad(m, 3, 0, where_str)
m2 = mb_pad(m, [1,2], 3, where_str)
m2 = mb_pad(m, 2, 'replicate', where_str)
m2 = mb_pad(m, 2, 'circular', where_str)
m2 = mb_pad(m, 2, 'symmetric', where_str)
m2 = mb_pad(m, [2,1], 'symmetric_repeat', where_str)
m2 = mb_pad(m, [-1,0], 'symmetric', where_str)
m2 = mb_pad(m, [0,-2], 'symmetric', where_str)
m2 = mb_pad(m, -1, '', where_str)
m2 = mb_pad(m, [2, -2], 3, where_str)
m


