function [f_tgt] = add_frame(f_src, f_tgt, offset, theta, fill_val)
%ADD_FRAME add a source frame to the target frame(s) given an offset and
%rotation angle theta
%   [frames] = add_frame(f_src, f_tgt, offset, theta)
%
% Inputs:
%      f_src - source frame
%
%      f_tgt - target frame
%
%      offset - off set in [x y] directions
%
%      theta - scalar in degrees
%
%
% Outputs:
%      frames - *Insert description of input variable here*
%
%
% Example:
%
% Notes: Offset and theta is relative to axe with an origin at the centre of
% each frame. Currently implemented so that the output frame size is the
% same as the target frame. Thus pixels defined in the target but not in
% the transformed source are inlcluded, but pixels in the transformed
% source not defined in the target are discarded
%
% Note the original size of the source frame should be the same as the
% original size of the target frame (the target frame may have subsequently
% been enlarged having had other src frames added)
%
% See also:
%
% Created: 24-Aug-2011
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

%set default to add to the back of the array
if nargin < 5
    fill_val = NaN;
end

%Rotate the source frame
if theta
    f_src = imrotate(f_src, theta, 'bicubic', 'crop');
end

%Now pad the src frame according to the offset
if offset(1) < 0
    f_src = padarray(f_src, [0 2*abs(offset(1))], fill_val, 'post');
elseif offset(1) > 0
    f_src = padarray(f_src, [0 2*offset(1)], fill_val, 'pre');
end
if offset(2) < 0
    f_src = padarray(f_src, [2*abs(offset(2)) 0], fill_val, 'post');
elseif offset(2) > 0
    f_src = padarray(f_src, [2*offset(2) 0], fill_val, 'pre');
end

%Post pad the target or source frame(s) as necessary - there should always
%be an even difference between dimensions of the src and target frame, so
%padding as applied to both sides of the respective frame
[rt ct nt] = size(f_tgt); 

%Get size of source following its transformation
[rs cs ns] = size(f_src);

if rt > rs
    f_src = padarray(f_src, [(rt - rs)/2 0], fill_val);
elseif rs > rt
    f_tgt = padarray(f_tgt, [(rs - rt)/2 0], fill_val);
end
if ct > cs
    f_src = padarray(f_src, [0 (ct - cs)/2], fill_val);
elseif cs > ct
    f_tgt = padarray(f_tgt, [0 (cs - ct)/2], fill_val);
end

% Match greylevels of corresponding pixels
tgt_ref = f_tgt(:,:,1);
inds = (~isnan(f_src) & ~isnan(tgt_ref));
src_pixels = f_src(inds);
tgt_pixels = tgt_ref(inds);
gain_offset = [src_pixels ones(size(src_pixels))] \ tgt_pixels;
f_src = f_src*gain_offset(1) + gain_offset(2);

%Finally we can add the source frame to the end of the target frame(s)
f_tgt = cat(3, f_tgt, f_src);





