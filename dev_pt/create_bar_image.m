function [prop_image] = create_bar_image(halfwidth, halfsize, theta_rad);

f_debug = (nargin == 0 && nargout == 0);
if f_debug
end

if ~exist('halfwidth','var'), halfwidth = 4; end
if ~exist('halfsize','var'), halfsize = 24; end
if ~exist('theta_rad','var'), theta_rad = 30 * pi/180; end

v = [-sin(theta_rad); cos(theta_rad)]; % normal to the line

fullsize = 2*halfsize+1;
samples = zeros(fullsize, fullsize);
pos_samples = zeros(fullsize, fullsize);
prop_image = zeros(fullsize, fullsize);

count = 0;
while (count < 10000)
    x = ((2*rand)-1) * halfsize;
    y = ((2*rand)-1) * halfsize;
    if (x*x+y*y > halfsize*halfsize)
        continue;
    end

    p = [x; y];
    d = p'*v;
    
    p = round(p);
    p = halfsize+1 + [p(1), -p(2)];
    samples(p(2), p(1)) = samples(p(2), p(1)) + 1;
    if (abs(d) < halfwidth)
        pos_samples(p(2), p(1)) = pos_samples(p(2), p(1)) + 1;
        count = count + 1;
    end
end
prop_image = pos_samples ./ samples;
prop_image(isnan(prop_image)) = 0;

if f_debug
    figure(1); clf; colormap(gray(256));
    imagesc(prop_image); axis('image','ij');
    clear prop_image;
end
