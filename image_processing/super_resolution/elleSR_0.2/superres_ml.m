function [hi_res_img] = superres_ml(W, Y, lambda1, lambda2, avg_img, num_itrs)
%
% function [hi_res_img,opts,flog,plog,sclog] = superres_ml(W,Y,La,Lb,avg_img,opts);
%
% W: W matrix (N*KM)
% Y: Y vector (KM*1)
% La: multiplicative photometric parameter.
% Lb: additive photometric parameter.
% avg_img: average image (biv*bih)
% opts: optional options vector to pass to scg.
%
if ~exist('num_itrs', 'var') || isempty(num_itrs)
    num_itrs = 20;
end

[ni_hi, nj_hi] = size(avg_img);

Y = Y-lambda2;
KM = numel(lambda1);
W = sparse(1:KM,1:KM,lambda1,KM,KM)*W;

x = avg_img(:);

options = optimoptions('fminunc',...
    'Algorithm','trust-region',...
    'Display', 'iter',...
    'SpecifyObjectiveGradient',true,...
    'HessianMultiplyFcn',@(H,Y)Y,...
    'MaxIterations', num_itrs);

tic;
obj_fun = @(x)likelihood(x,...
    W, Y);
X = fminunc(obj_fun, x,...
    options);
toc;
hi_res_img = reshape(X,ni_hi,nj_hi);


% GRADTEST = 0; % Set this to 1 in order to run local_gradcheck below.
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
% 
% biv = size(avg_img,1);
% bih = size(avg_img,2);
% 
% if debug
%     debug_im_sz = [biv, bih];
% else
%     debug_im_sz = 0;
% end
% 
% Y = Y-lambda2;
% KM = numel(lambda1);
% 
% %W = mex_amub(W,sparse(1:KM,1:KM,La));
% % The line above would be "W = W*sparse(1:KM,1:KM,La);", but matlab doesn't
% % handle memory very well in this case, so it generates a completely
% % unnecessary "out of memory" error. I got mex_amub from the file exchange,
% % though now "ssmult" seems to be preferred.
% W = W*sparse(1:KM,1:KM,lambda1,KM,KM);
% 
% x = avg_img(:)';
% 
% if (GRADTEST)
%     [D,G] = local_gradcheck(x,W,Y);
%     hi_res_img = D; opts = G; % overwrite usual outputs when testing grad.
%     return;
% end
% 
% if (opts(1)>=0), fprintf('hmem_ml iterating (%i iters):\n',opts(14)); end;
% if (nargout <4)
%     [X,opts,flog] = scg(@local_func,x,opts,@local_grad,W,Y,debug_im_sz);
% else
%     [X,opts,flog,plog,sclog] = scg(@local_func,x,opts,@local_grad,W,Y,debug_im_sz);
% end
% hi_res_img = reshape(X,biv,bih);
% if (opts(1)>=0), fprintf('Done!\n'); end;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [L, Ldx, H] = likelihood(x, W, Y)

r = Y - W*x; %W is alreday pre-multiplied by lambda1
L = sum(r.^2);

if nargout > 1
   Ldx = -2*(W'*r); 
end
if nargout > 2
    H = [];
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e = local_func(x,W,Y,debug_im_sz)
if debug_im_sz(1)
    figure; imgray(reshape(x,debug_im_sz));
end
e = Y'-x*W;
e = sum(e.^2);
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function g = local_grad(x,W,Y,debug) %#ok
r = Y'-x*W;
g = W*r';
g = -2.*g';
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