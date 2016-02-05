function complex_out = rgb2complex(c_rgb, phase_lims, mag_scale, max_size, bg_val)
% Map a complex number to the rgb colorspace using a transform into the hsv
% colorspace, such that complex phase is represented as hue, saturation is
% 1, and magnitude is represented as brightness

if (nargin==0 && nargout==0), test(); return; end

%Get size of image;
[row col dims] = size(c_rgb);

if dims ~= 3
    error('Input image must by RGB');
end

%pre-allocate complex_out
complex_out = zeros(row, col);

%map phase from limits to [0, 1)
if ~exist('phase_lims', 'var') || isempty(phase_lims), phase_lims = [-pi pi]; end
phase_scale = phase_lims(2) - phase_lims(1);
phase_offset = phase_lims(1);

%map mag to [0, 1)
if ~exist('mag_scale', 'var') || isempty(mag_scale), mag_scale =  1; end

%compute number of parts we need to break image into
if ~exist('max_size', 'var') || isempty(max_size), max_size = 512; end
r_parts = ceil(row / max_size);
c_parts = ceil(col / max_size);

if ~exist('bg_val', 'var') || isempty(bg_val), bg_val = 0; end;

%Go through each segment computing hsv2rgb
for rp = 1:r_parts
    for cp = 1:c_parts
        sr = 1 + (rp-1)*max_size;
        er = min(rp*max_size, row);
        sc = 1 + (cp-1)*max_size;
        ec = min(cp*max_size, col);
        
        %Generate HSV representation for this part
        c_hsv = rgb2hsv(c_rgb(sr:er,sc:ec,:));
        c_phase = c_hsv(:,:,1)*phase_scale + phase_offset;
        if bg_val
            c_mag = mag_scale * (c_hsv(:,:,2) - 1 + bg_val) / bg_val;
        else
            c_mag = mag_scale * c_hsv(:,:,3);
        end
        
        %Convert to rgb
        complex_out(sr:er,sc:ec) = c_mag .* exp(1.0i * c_phase);
    end
end


%% Test function
function test()
clc;

sz = 75;
r = linspace(-1,1, sz);
[xx,yy] = meshgrid(r,r);
rr = sqrt(xx.^2 + yy.^2);

xx(rr > 1.0) = 0;
yy(rr > 1.0) = 0;

rgb = complex2rgb(complex(xx,yy), [],[],[], 1);
image(rgb);
axis('image','ij');

imwrite(rgb, 'U:\projects\nailfold\synthesis\showcase\create_flowmap\flow_legend.png');


