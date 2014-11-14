function c_rgb = complex2rgb(complex_in, phase_lims, mag_lim, max_size, bg_val)
% Map a complex number to the rgb colorspace using a transform into the hsv
% colorspace, such that complex phase is represented as hue, saturation is
% 1, and magnitude is represented as brightness

if (nargin==0 && nargout==0), test(); return; end

%Get size of image;
[row col] = size(complex_in);

%pre-allocate c_rgb
c_rgb = zeros(row, col, 3);

%Take phase and magnitude of complex input
c_phase = angle(complex_in);
c_mag = abs(complex_in);
clear complex_in;

%map phase from limits to [0, 1)
if ~exist('phase_lims', 'var') || isempty(phase_lims), phase_lims = [-pi pi]; end
c_phase = (c_phase - phase_lims(1)) / (phase_lims(2) - phase_lims(1));
c_phase(c_phase < 0) = 0;
c_phase(c_phase > 1) = 1;
c_phase(isnan(c_phase)) = 0;

%map mag to [0, 1)
if ~exist('mag_lim', 'var') || isempty(mag_lim), mag_lim =  max(c_mag(:)); end
c_mag = c_mag ./ mag_lim;
c_mag(c_mag > 1) = 1;
c_mag(isnan(c_mag)) = 0;

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
        c_hsv        = c_phase(sr:er,sc:ec,:);
        c_hsv(:,:,2) = (1-bg_val) + bg_val*c_mag(sr:er,sc:ec,:);
        c_hsv(:,:,3) = bg_val + (1-bg_val)*c_mag(sr:er,sc:ec,:);
        
        %Convert to rgb
        c_rgb(sr:er,sc:ec,:) = hsv2rgb(c_hsv);
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


