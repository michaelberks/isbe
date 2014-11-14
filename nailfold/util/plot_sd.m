function plot_sd(x,y,y_plus,y_minus,color,alpha)

f_debug = (nargin==0) && (nargout==0);
if f_debug
	x = linspace(0,2*pi,100);
	y = sin(x);
	y_plus = y + 0.2;
	y_minus = y - 0.2;
end

if nargin<5, color = 'b'; end;
if nargin<6, alpha = 0.1; end; % strength of background shading

hold on;
h = patch(x([1:end,end:-1:1]),[y_plus,y_minus(end:-1:1)],color);
col = get(h,'facecolor');
lightcolor = (1-alpha) + (alpha*col);
set(h,'facecolor',lightcolor,'edgecolor',lightcolor);
plot(x,y,color);

if f_debug
	clear;
end
