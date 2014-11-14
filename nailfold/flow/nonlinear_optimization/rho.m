function err_out = rho(err_in, fname, varargin)
% Apply robust error function to error vector in order to downweight
% outliers

if (nargin==0 && nargout==0), test(); return; end

if ~exist('fname','var') || isempty(fname), fname = ''; end

switch fname
    case {'lorentzian'},
        if length(varargin)>0, nu = varargin{1}; end
        if ~exist('nu','var'), nu = 1; end

        err_out = log(1 + 0.5*(err_in/nu).^2);
        
        if any(isnan(err_out(:))) || ...
           any(isinf(err_out(:))) || ...
           any(err_out(:) < 0) || ...
           any(~isreal(err_out(:)))
            keyboard;
        end
        
    case {'charbonnier'},
        if length(varargin)>0, epsilon = varargin{1}; end
        if ~exist('epsilon','var'), epsilon = 0.001; end

        if length(varargin)>1, a = varargin{2}; end
        if ~exist('a','var'), a = 0.45; end
        
        err_out = (err_in.^2 + epsilon).^a - epsilon^a;
                
    case {'geman-mcclure'},
        if length(varargin)>0, mu = varargin{1}; end
        if ~exist('mu','var'), mu = 1; end
        
        e2 = err_in.^2;
        err_out = e2 ./ (e2 + mu*mu);

    case {'quadratic'}
        % Return squared error
        err_out = err_in.^2;
        
    case {'linear'}
        % Return error unchanged
        err_out = err_in;
        
    otherwise,
        error(['Unknown error function: ', fname]);
end


%% Test function
function test()
clc;

xmax = 10;
x = linspace(-xmax, xmax, 101);

legs = {'quadratic','lorentzian','charbonnier','geman-mcclure'};

y = [];
for i = 1:length(legs)
    y(i,:) = rho(x, legs{i});
end

% Close any existing figures from this function, without the error if none
% exist.
try close('fig_rho'); catch end

figure; clf; hold on;
set(gcf,'name','fig_rho');
    plot(x, y);
    legend(legs,'location','north');
    axis([-xmax,xmax, 0,xmax]);
    