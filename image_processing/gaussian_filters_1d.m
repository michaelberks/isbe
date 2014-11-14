function [g,dg,ddg,x] = gaussian_filters_1d(varargin)
% aux function to make the 1d gaussian derivative filters

f_debug = (nargin==0 && nargout==0);
if f_debug, test_script(); return; end

[g,dg,ddg,x] = func(varargin{:});


%% The function
function [g,dg,ddg,x] = func(sigma, width)
if ~exist('sigma','var'), sigma = 1; end
if ~exist('width','var'), width = round(5*sigma); end

sigmasq = sigma^2;

x	= (-width:width);
k	= (2*pi)^-0.5; %Normalising factor without sigma
%k	= (2*pi*sigmasq)^-0.5; %Normalising factor with sigma

g	= k*exp(-0.5* (x.*x)/sigmasq);
dg	= -x/sigmasq .* g;
ddg	= (-1/sigmasq * g) - (x/sigmasq .* dg);

% Shift thefirst and second derivatives to have zero sum
% (will make little difference for sufficient filter width)
dg = dg - mean(dg);
ddg = ddg - mean(ddg);

%Don't do this!!
% g	= g / sum(abs(g));
% dg	= dg / sum(abs(dg));
% ddg	= ddg / sum(abs(ddg));


%% Test script
function test_script()

sigma = 4;
width = 5*sigma;

[g,dg,ddg,x] = func(sigma, width);

figure(1); clf; hold on;
plot(x,g,'b-');
plot(x,dg,'r-');
plot(x,ddg,'g-');
clear;

return

for i = 1:6
    switch i
        case 1, im = g'*dg; filename = 'Gx.png';
        case 2, im = dg'*g; filename = 'Gy.png';
        case 3, im = g'*ddg; filename = 'Gxx.png';
        case 4, im = ddg'*g; filename = 'Gyy.png';
        case 5, im = g'*ddg-ddg'*g; filename = 'Gxx-Gyy.png';
        case 6, im = dg'*dg; filename = 'Gxy.png';
    end

% 		im = im-min(im(:));
% 		im = im/max(im(:));
    im = uint8(127*(1+im/max(abs(im(:)))));

    figure(1); clf; colormap(gray(256));
    image(im);

    axis('image','off');
    imwrite(im,['S:\projects\mammography\matlab\papers\2011cvpr\figs\',filename]);
end
clear;
