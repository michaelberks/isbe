function [hi_res_img] = superres_huber(W, Y, lambda1, lambda2, avg_im, alpha, beta, num_itrs)
%
% function [im,opts,flog,plog,sclog] = superres_huber(W,Y,La,Lb,avg_im,alpha,beta,opts);
%
% W: W matrix (N*KM)
% Y: Y vector (KM*1)
% La: multiplicative photometric parameter.
% Lb: additive photometric parameter.
% avg_im: average image (biv*bih)
% alpha: Huber alpha parameter.
% beta: (beta) -- ratio of Huber prior strength to image noise precision.
% opts: optional options vector to pass to scg.
%
% if (nargin<6)
%     fprintf('Using default options vector.\n');
%     opts = zeros(1,18);
%     opts(1) = 1;
%     opts(2) = 1e-10;
%     opts(3) = 1e-10;
%     opts(14) = 40;
% end
% if nargin < 7
%     debug = false;
% end
if ~exist('num_itrs', 'var')
    num_itrs = 20;
end

[ni_hi, nj_hi] = size(avg_im);

Y = Y-lambda2;
KM = numel(lambda1);
W = sparse(1:KM,1:KM,lambda1,KM,KM)*W;

x = avg_im(:);

% options = optimoptions('fminunc',...
%     'Algorithm','trust-region',...
%     'Display', 'iter',...
%     'SpecifyObjectiveGradient',true,...
%     'HessianMultiplyFcn',@(H,Y)Y,...
%     'MaxIterations', num_itrs);
% 
% tic;
% obj_fun = @(x)likelihood(x,...
%     W, Y, alpha, beta, ni_hi, nj_hi);
% X = fminunc(obj_fun, x,...
%     options);
% toc;

opts = zeros(1,18);
opts(1) = 1;
opts(2) = 1e-10;
opts(3) = 1e-10;
opts(14) = num_itrs;
[X] = scg(@local_func, x, opts, @local_grad, W, Y, alpha, beta, ni_hi, nj_hi);

hi_res_img = reshape(X,ni_hi,nj_hi);

% tic;
% if (opts(1)>=0), fprintf('hmem_ml iterating (%i iters):\n',opts(14)); end;
% if (nargout <4)
%     [X,opts,flog] = scg(@local_func,x,opts,@local_grad,W,Y,alpha,beta,ni_hi,nj_hi,debug);
% else
%     [X,opts,flog,plog,sclog] = scg(@local_func,x,opts,@local_grad,W,Y,alpha,beta,ni_hi,nj_hi,debug);
% end
% toc;
% im = reshape(X,ni_hi,nj_hi);
% 
% 
% if (opts(1)>=0), fprintf('Done!\n'); end;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [L Ldx H] = likelihood(x,W,Y,alpha,beta,ni_hi,nj_hi)

r = Y - W*x; %W is alreday pre-multiplied by lambda1
Dx = abs(diff_x(x, ni_hi, nj_hi));

L = sum(r.^2) + beta*huber_prior(Dx,alpha);

if nargout > 1
    Ldx = -2*(W'*r) + beta*huber_prior_grad(Dx, alpha, ni_hi, nj_hi);
end
if nargout > 2
    H = [];
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L = local_func(x,W,Y,alpha,beta,ni_hi,nj_hi)

r = Y - W*x;
Dx = abs(diff_x(x, ni_hi, nj_hi));
L = sum(r.^2) +  beta*huber_prior(Dx,alpha);
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Ldx = local_grad(x,W,Y,alpha,beta,ni_hi,nj_hi)
r = Y - W*x;
Dx = abs(diff_x(x, ni_hi, nj_hi));
Ldx = -2*(W'*r) + beta*huber_prior_grad(Dx, alpha, ni_hi, nj_hi);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D,G] = local_gradcheck(x,varargin)
% This function checks the gradient numerically, for a random subsection of
% the image pixels (since it takes rather too long to do all of them as in
% scg's automatic gradient test option, selected by setting opts(9)=1).

ep = 1e-4; % numerical grdient step size.
numsamp = 30; % approximate number of pixels I'd like to check.

x = x + reshape(gsamp(0,0.15^2,numel(x)),size(x)); % Add in a bit of randome noise.
G = local_grad(x,varargin{:});
L = numel(x);
if (L~=numel(G)), error('lengths of parameter vector and gradient vector are inconsistent.'); end;
step = zeros(1,L);
D = step;
t = numsamp/L;
for i = 1:L
    if (mod(i,10000)<1), fprintf('Element %i of %i.\n',i,L); end;
    if(rand<t)
        step(i) = ep;
        l1 = local_func(x+step,varargin{:});
        l2 = local_func(x-step,varargin{:});
        D(i) = (0.5/ep)*(l1-l2);
        step(i) = 0;
    else
        D(i) = NaN;
    end
end
fprintf('Done all elements!\n');
m = logical(~isnan(D));
G = G(m);
D = D(m);
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%