function [hi_res_img, lambda_new, theta_new, alpha_new, beta_new] = ...
    super_res_simultaneous(Y, theta, lambda, gamma, ni_lo, nj_lo, ni_hi, nj_hi, varargin)
%SUPER_RES_SIMULTANEOUS *Insert a one line summary here*
%   [] = super_res_simultaneous(varargin)
%
% SUPER_RES_SIMULTANEOUS uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%   avg_im: average image (ni_hi*nj_hi)
%   W: W matrix (N*KM)
%   Y: Y vector (KM*1)
%   La: multiplicative photometric parameter.
%   Lb: additive photometric parameter.
%
% Optional Arguments:
%   alpha: Huber args.alpha parameter.
%   beta: (args.beta) -- ratio of Huber prior strength to image noise precision.
%   opts: optional options vector to pass to scg.
%   do_grad_test:
%   debug:
%
% Outputs:
%   hi_res_img: Super-resolution image
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 08-Mar-2017
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin, 0, ...
    'alpha', 0.01,...
	'beta', 0.1,...
    'theta_offset', [],...
    'theta_scaling', [], ...
    'theta_gt', [],...
    'lambda_offset', [],...
    'lambda_scaling', [], ...
    'n_itr_all', 2,...
	'n_itr_alpha', 10,...
    'n_itr_x', 20,...
	'do_grad_test', false,...
    'validation_pct', 0.05,...
    'x_thresh', 0.3,...
    'lambda_thresh', 0.3,...
    'theta_thresh', 0.3,...
    'alpha_thresh', 10e-4,...
    'beta_thresh', 10e-4,...
    'debug', false);
clear varargin;

%Get K and size of images;
K = size(theta, 2);
if length(ni_lo) == 1
    ni_lo = ni_lo*ones(K,1);
end
if length(nj_lo) == 1
    nj_lo = nj_lo*ones(K,1);
end
if length(gamma) == 1
    gamma = gamma*ones(K,1);
end

if isempty(args.theta_offset)
    args.theta_offset = mean(theta, 2);
end
if isempty(args.theta_scaling)
    args.theta_offset = 0.35 / std(theta, 1, 2);
end

if isempty(args.lambda_offset)
    args.lambda_offset = mean(lambda, 2);
end
if isempty(args.lambda_scaling)
    args.lambda_offset = 0.35 / std(lambda, 1, 2);
end

%1) Initialise x0, L_0, theta0, alpha0, beta0
alpha0 = args.alpha;
beta0 = args.beta;
[theta0, lambda0, W, avg_img] = ...
    initialise_theta_lambda(ni_hi, nj_hi, ni_lo, nj_lo, Y, theta, lambda, gamma, ceil(K/4),...
    args.theta_offset, args.theta_scaling, args.lambda_offset, args.lambda_scaling);

[lambda1, lambda2] = expand_lambda(lambda0, ni_lo, nj_lo);
[x0] = initialise_x(W, Y, lambda1, lambda2, avg_img, ceil(K/4));

%Now enter the main outer optimisation loop
converged = false;
itr = 0;
while itr < args.n_itr_all && ~converged
    
%     %2a) Sample set of validation pixels
%     [val_px] = sample_validation_pixels(ni_lo, nj_lo, args.validation_pct);
%     
%     %2b) Update alpha and beta
%     [alpha_new, beta_new] = ...
%        update_alpha_beta(x0, theta0, lambda0, alpha0, beta0,...
%        Y, W, val_px, gamma, ni_hi, nj_hi, ni_lo, nj_lo, args.n_itr_alpha);
    alpha_new = alpha0;
    beta_new = beta0;
    
    %2c) Update x, L and theta
    [x_new, theta_new, lambda_new, W] = ...
        update_x_theta_lambda(x0, Y, theta0, lambda0, alpha_new, beta_new,...
        gamma, ni_hi, nj_hi, ni_lo, nj_lo, args.n_itr_x,...
        args.theta_offset, args.theta_scaling, args.lambda_offset, args.lambda_scaling);
    
    %3) Check convergence
    converged = check_converged(....
        x_new, lambda_new, theta_new, alpha_new, beta_new,...
        x0, lambda0, theta0, alpha0, beta0,...
        args.x_thresh, args.lambda_thresh, args.theta_thresh, args.alpha_thresh, args.beta_thresh);
    
    %4)Move new to old for start of next iteration
    if (~converged)
        x0 = x_new;
        lambda0 = lambda_new;
        theta0 = theta_new;
        alpha0 = alpha_new;
        beta0 = beta_new;
        itr = itr + 1;
    end
    
end

%Reshape the final high resolution and return final values for the
%photometric, geometric and noise paramaters
hi_res_img = reshape(x_new, ni_hi, nj_hi);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [theta0, lambda0, W, avg_img] = ...
    initialise_theta_lambda(ni_hi, nj_hi, ni_lo, nj_lo, Y, theta, lambda, gamma, num_itr,...
    theta_offset, theta_scaling, lambda_offset, lambda_scaling)

%Compute average image from initial data
[avg_img] = getAvim(ni_hi, nj_hi, ni_lo, nj_lo, Y, theta, gamma);

%Run a few optimisation steps to initialize L and theta using avg_img in
%place of x

options_theta = optimoptions('fminunc',...
    'Algorithm','trust-region',...quasi-newton
    'Display', 'iter',...
    'SpecifyObjectiveGradient',true,...false
    'HessianMultiplyFcn',@(H,Y)Y,...[]
    'MaxIterations', num_itr);

options_lambda = optimoptions('fminunc',...
    'Algorithm','trust-region',...
    'Display', 'iter',...
    'SpecifyObjectiveGradient',true,...
    'HessianMultiplyFcn',@(H,Y)Y,...
    'MaxIterations', num_itr);

[lambda1, lambda2] = expand_lambda(lambda, ni_lo, nj_lo);
x = avg_img(:);

%objective function, optimising (stack of photometric L, and geometric theta registration params) 
obj_fun = @(theta)likelihood_theta(theta,...
    lambda1, lambda2, x, Y, gamma,... %fixed variables in the objective function
    ni_lo, nj_lo, ni_hi, nj_hi,...
    theta_offset, theta_scaling); %auxilliary variables required to compute objective
theta = scale_theta(theta, theta_offset, theta_scaling);
theta0 = fminunc(obj_fun, theta(:),...
    options_theta); 
theta0 = rescale_theta(reshape(theta0, size(theta)), theta_offset, theta_scaling);

%Make new W using the new theta, clearing the old from memory first
W = makeW(ni_hi, nj_hi, ni_lo, nj_lo, theta0, gamma);

lambda = scale_lambda(lambda, lambda_offset, lambda_scaling);
obj_fun = @(lambda)likelihood_lambda(lambda,...
    x, W, Y,... %fixed variables in the objective function
    ni_lo, nj_lo,...
    lambda_offset, lambda_scaling); %auxilliary variables required to compute objective
lambda0 = fminunc(obj_fun, lambda(:),...
    options_lambda);
lambda0 = rescale_lambda(reshape(lambda0, size(lambda)), lambda_offset, lambda_scaling);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x0] = initialise_x(W, Y, lambda1, lambda2, avg_img, num_itr)
%Apply ML optimisation using estimated registration parameters
%and the average image
[x0] = superres_ml(W, Y, lambda1, lambda2, avg_img, num_itr);
x0 = x0(:);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_new, theta_new, lambda_new, W_new] = update_x_theta_lambda(x, Y, theta, lambda, alpha, beta, gamma,...
    ni_hi, nj_hi, ni_lo, nj_lo, num_itr, theta_offset, theta_scaling, lambda_offset, lambda_scaling)

options_theta = optimoptions('fminunc',...
    'Algorithm','trust-region',...quasi-newton
    'Display', 'iter',...
    'SpecifyObjectiveGradient',true,...false
    'HessianMultiplyFcn',@(H,Y)Y,...[]
    'MaxIterations', num_itr);

options_x_lambda = optimoptions('fminunc',...
    'Algorithm','trust-region',...
    'Display', 'iter',...
    'SpecifyObjectiveGradient',true,...
    'HessianMultiplyFcn',@(H,Y)Y,...
    'MaxIterations', num_itr);

[lambda1, lambda2] = expand_lambda(lambda, ni_lo, nj_lo);

%objective function, optimising stack of geometric registration params (theta) 
obj_fun = @(theta)likelihood_theta(theta(:),...
    lambda1, lambda2, x, Y, gamma,... %fixed variables in the objective function
    ni_lo, nj_lo, ni_hi, nj_hi,...
    theta_offset, theta_scaling); %auxilliary variables required to compute objective
theta = scale_theta(theta, theta_offset, theta_scaling);
theta_new = fminunc(obj_fun, theta(:),...
    options_theta); 
theta_new = rescale_theta(reshape(theta_new, size(theta)), theta_offset, theta_scaling);

%Make new W using the new theta, clearing the old from memory first
W_new = makeW(ni_hi, nj_hi, ni_lo, nj_lo, theta_new, gamma);

%objective function, optimising stack of hi-res image and photometric params (x, lambda) 
obj_fun = @(x_lambda)likelihood_x_lambda(x_lambda,...
    W_new, Y, alpha, beta,... %fixed variables in the objective function
    ni_hi, nj_hi, ni_lo, nj_lo,...
    lambda_offset, lambda_scaling); %auxilliary variables required to compute objective
lambda = scale_lambda(lambda, lambda_offset, lambda_scaling);
x_lambda = fminunc(obj_fun, [x; lambda(:)],...
    options_x_lambda);

N = ni_hi * nj_hi;
x_new = x_lambda(1:N);
lambda_new = rescale_lambda(reshape(x_lambda(N+1:end), size(lambda)), lambda_offset, lambda_scaling);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [alpha1, beta1] = update_alpha_beta(x, theta, lambda, alpha, beta, Y, W, val_px, gamma, ni_hi, nj_hi, ni_lo, nj_lo, num_itr)

%Apply photometric params
[lambda1, lambda2] = expand_lambda(lambda, ni_lo, nj_lo);
Y = Y-lambda2;
KM = numel(lambda1);
W = sparse(1:KM,1:KM,lambda1,KM,KM)*W;

%Discard pixels from W and Y
W(val_px(:),:) = [];
y_val0 = Y(val_px(:));
Y(val_px(:)) = [];

options = optimoptions('fminunc',...
    'Display', 'iter',...
    'Algorithm','quasi-newton',...
    'SpecifyObjectiveGradient',false,...
    'MaxIterations', num_itr);

obj_fun = @(alpha_beta)compute_alpha_beta_err(alpha_beta,...
    x, theta, lambda , Y, W, y_val0, val_px, gamma, ni_hi, nj_hi, ni_lo, nj_lo, 3);
alpha_beta1 = fminunc(obj_fun, [log(alpha); log(beta)],...
    options);

alpha1 = exp(alpha_beta1(1));
beta1 = exp(alpha_beta1(2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y_err] = compute_alpha_beta_err(alpha_beta, x, theta, lambda, Y, W, y_val0, val_px, gamma, ni_hi, nj_hi, ni_lo, nj_lo, num_itr)

%Update x at the new alpha and beta
options = optimoptions('fminunc',...
    'Algorithm','trust-region',...
    'Display', 'off',...
    'SpecifyObjectiveGradient',true,...
    'HessianMultiplyFcn',@(H,Y)Y,...
    'MaxIterations', num_itr);

obj_fun = @(x)likelihood_x(x,...
    W, Y, exp(alpha_beta(1)), exp(alpha_beta(2)), ni_hi, nj_hi);
x1 = fminunc(obj_fun, x(:),...
    options);
x1 = reshape(x1, ni_hi, nj_hi);

%Re-estimate each low res image and sample the desired pixels...
y_val1 = zeros(size(y_val0));
curr_val_pt = 0;
curr_k_pt = 0;
for i_k = 1:length(ni_lo)
    n_k_px = ni_lo(i_k)*nj_lo(i_k);
    k_idx = curr_k_pt + (1:n_k_px);
    curr_k_pt = curr_k_pt + n_k_px;
    val_px_k = val_px(k_idx);
    
    n_val_px = sum(val_px_k);
    val_idx = curr_val_pt + (1:n_val_px);
    curr_val_pt = curr_val_pt + n_val_px;
    
    [lo_res_img] = makeLR(x1, ni_lo(i_k), nj_lo(i_k), theta(:,i_k), lambda(:,i_k), 0, gamma(i_k), val_px_k);
    y_val1(val_idx) = lo_res_img(val_px_k);
end

y_err = mean(abs(y_val0 - y_val1));

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [L Ldx H] = likelihood_x(x,W,Y,alpha,beta,ni_hi,nj_hi)

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
function [L, Ldxl, H] = likelihood_x_lambda(x_lambda, W, Y, alpha, beta, ni_hi, nj_hi, ni_lo, nj_lo, lambda_offset, lambda_scaling)

N = ni_hi*nj_hi;
K = length(ni_lo);

x = x_lambda(1:N);
lambda = rescale_lambda(reshape(x_lambda(N+1:end), 2, K), lambda_offset, lambda_scaling);
[lambda1, lambda2] = expand_lambda(lambda, ni_lo, nj_lo);

r = Y - lambda1.*(W*x) - lambda2;
Dx = abs(diff_x(x, ni_hi, nj_hi));

L = sum(r.^2) + beta*huber_prior(Dx,alpha);

if nargout > 1
    
    Ldx = -2*(W'*r) + beta*huber_prior_grad(Dx, alpha, ni_hi, nj_hi);
    
    Ldl = zeros(2, K);
    curr_pt = 0;
    for i_k = 1:K              
        n_pts = ni_lo(i_k)*nj_lo(i_k);
        idx = curr_pt + (1:n_pts);
        curr_pt = curr_pt + n_pts;
        r_k = r(idx);
        W_k = W(idx,:);
        
        Ldl(1, i_k) = -2*x'*W_k'*r_k;
        Ldl(2, i_k) = sum(r_k);          
    end
    Ldl = scale_lambda(Ldl, lambda_offset, lambda_scaling);
    Ldxl = [Ldx; Ldl(:)];
end

if nargout > 2
    H = [];
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [L, Ldl, H] = likelihood_lambda(lambda, x, W, Y, ni_lo, nj_lo, lambda_offset, lambda_scaling)

K = length(ni_lo);
lambda = rescale_lambda(reshape(lambda, 2, K), lambda_offset, lambda_scaling);
[lambda1, lambda2] = expand_lambda(lambda, ni_lo, nj_lo);
r = Y - lambda1.*(W*x) - lambda2;
L = sum(r.^2);

if nargout > 1
    
    Ldl = zeros(2, K);
    curr_pt = 0;
    for i_k = 1:K              
        n_pts = ni_lo(i_k)*nj_lo(i_k);
        idx = curr_pt + (1:n_pts);
        curr_pt = curr_pt + n_pts;
        r_k = r(idx);
        W_k = W(idx,:);
        
        Ldl(1, i_k) = -2*x'*W_k'*r_k;
        Ldl(2, i_k) = sum(r_k);          
    end
    Ldl = scale_lambda(Ldl, lambda_offset, lambda_scaling);
    Ldl = Ldl(:);
end
if nargout > 2
    H = [];
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [L, Ldtheta, H] = likelihood_theta(theta, lambda1, lambda2, x, Y, gamma,...
    ni_lo, nj_lo, ni_hi, nj_hi, theta_offset, theta_scaling)

K = length(ni_lo);
theta = rescale_theta(reshape(theta, 8, K), theta_offset, theta_scaling);
W = makeW(ni_hi, nj_hi, ni_lo, nj_lo, theta, gamma);
Y = Y - lambda2;
r = Y - lambda1.*(W*x);
L = sum(r.^2);

if nargout > 1
    Ldtheta = zeros(8, K);
    curr_pt = 0;
    for i_k = 1:K 
        n_pts = ni_lo(i_k)*nj_lo(i_k);
        idx = curr_pt + (1:n_pts);
        curr_pt = curr_pt + n_pts;
        Y_k = Y(idx);
        lambda1_k = lambda1(idx);
        
        %Sub-sample each of the above at a set of validations pixels
        [val_px] = sample_validation_pixels(ni_lo(i_k), nj_lo(i_k), 0.5);           
        delta = 0.001;
        for i_th = 1:8
            theta_p = theta;
            theta_n = theta;
            theta_p(i_th, i_k) = theta_p(i_th, i_k) + delta;
            theta_n(i_th, i_k) = theta_n(i_th, i_k) - delta;
            W_kp = update_Wk(ni_hi, nj_hi, ni_lo(i_k), nj_lo(i_k), theta_p(:,i_k), gamma(i_k), val_px);
            W_kn = update_Wk(ni_hi, nj_hi, ni_lo(i_k), nj_lo(i_k), theta_n(:,i_k), gamma(i_k), val_px);

            r_kp = Y_k - lambda1_k.*(W_kp*x);
            L_kp = sum(r_kp(val_px).^2);
            r_kn = Y_k - lambda1_k.*(W_kn*x);
            L_kn = sum(r_kn(val_px).^2);
            
            Ldtheta(i_th, i_k) = (L_kp - L_kn) / (2*delta);

        end
            
    end
    Ldtheta = scale_theta(Ldtheta, theta_offset, theta_scaling);
    Ldtheta = Ldtheta(:);
end
if nargout > 2
    H = [];
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [val_px] = sample_validation_pixels(ni_lo, nj_lo, validation_pct)
KM = sum(ni_lo.*nj_lo);
val_px = false(KM,1);
K = length(ni_lo);
curr_pt = 0;
for i_k = 1:K
    n_pts = ni_lo(i_k) * nj_lo(i_k); 
    n_samples = floor(n_pts * validation_pct);
    
    [potential_cols, potential_rows] = meshgrid(5:nj_lo(i_k)-4, 5:ni_lo(i_k)-4);
    potential_idx = sub2ind([ni_lo(i_k) nj_lo(i_k)], potential_rows(:), potential_cols(:));
    r_samples = randperm(length(potential_idx), n_samples);
    
    sample_idx = curr_pt + potential_idx(r_samples);
    val_px(sample_idx) = 1;
    curr_pt = curr_pt + n_pts;
end
    
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [theta_t] = scale_theta(theta, theta_offset, theta_scaling)
theta_t = bsxfun(@times, bsxfun(@minus, theta, theta_offset), theta_scaling);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [theta] = rescale_theta(theta_t, theta_offset, theta_scaling)
theta = bsxfun(@plus, bsxfun(@rdivide, theta_t, theta_scaling), theta_offset);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lambda_t] = scale_lambda(lambda, lambda_offset,lambda_scaling)
lambda_t = bsxfun(@times, bsxfun(@minus, lambda, lambda_offset), lambda_scaling);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lambda] = rescale_lambda(lambda_t, lambda_offset,lambda_scaling)
lambda = bsxfun(@plus, bsxfun(@rdivide, lambda_t, lambda_scaling), lambda_offset);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lambda1, lambda2] = expand_lambda(lambda, ni_lo, nj_lo)
KM = sum(ni_lo .* nj_lo);
lambda1 = ones(KM,1);
lambda2 = ones(KM,1);
curr_pt = 0;
for i_k = 1:length(ni_lo)
    n_pts = ni_lo(i_k)*nj_lo(i_k);
    idx = curr_pt + (1:n_pts);
    curr_pt = curr_pt + n_pts;
    lambda1(idx) = lambda(1,i_k);
    lambda2(idx) = lambda(2,i_k);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lambda] = collapse_lambda(lambda1, lambda2, ni_lo, nj_lo)%#ok
K = length(ni_lo);

lambda = ones(2, K);
curr_pt = 1;
for i_k = 1:length(ni_lo)
    
    lambda(1,i_k) = lambda1(curr_pt);
    lambda(2,i_k) = lambda2(curr_pt);
    
    n_pts = ni_lo(i_k)*nj_lo(i_k);
    curr_pt = curr_pt + n_pts;

end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [converged] = check_converged(....
        x_new, lambda_new, theta_new, alpha_new, beta_new,...
        x, lambda, theta, alpha, beta,...
        x_thresh, lambda_thresh, theta_thresh, alpha_thresh, beta_thresh)
converged = ...
    all(abs(x_new(:)      - x(:))       < x_thresh) && ...
    all(abs(lambda_new(:) - lambda(:))  < lambda_thresh) && ...
    all(abs(theta_new(:)  - theta(:))   < theta_thresh) && ...
    all(abs(alpha_new(:)  - alpha(:))   < alpha_thresh) && ...
    all(abs(beta_new(:)   - beta(:))    < beta_thresh);
return

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [L, Ldtl, H] = likelihood_theta_lambda(theta_lambda, x, W, Y, gamma, ni_lo, nj_lo, ni_hi, nj_hi) %#ok
% 
% K = length(ni_lo);
% KM = size(Y,1);
% theta = reshape(theta_lambda(1:8*K), 8, K);
% lambda = reshape(theta_lambda(8*K+1:end), 2, K);
% 
% lambda1 = ones(KM,1);
% lambda2 = ones(KM,1);
% curr_pt = 0;
% for i_k = 1:K
%     n_pts = ni_lo(i_k)*nj_lo(i_k);
%     idx = curr_pt + (1:n_pts);
%     curr_pt = curr_pt + n_pts;
%     lambda1(idx) = lambda(1,i_k);
%     lambda2(idx) = lambda(2,i_k);
% end
% 
% r = Y - lambda1.*(W*x) - lambda2;
% L = sum(r.^2);
% 
% if nargout > 1
%     
%     Ldtheta = zeros(8, K);
%     Ldl1 = zeros(K,1);
%     Ldl2 = zeros(K,1);
%     curr_pt = 0;
%     for i_k = 1:K              
%         n_pts = ni_lo(i_k)*nj_lo(i_k);
%         idx = curr_pt + (1:n_pts);
%         curr_pt = curr_pt + n_pts;
%         r_k = r(idx);
%         W_k = W(idx,:);
%         L_k = sum(r_k.^2);
%         Y_k = Y(idx);
%         lambda1_k = lambda1(idx);
%         lambda2_k = lambda2(idx);
%         
%         Ldl1(i_k) = -2*x'*W_k'*r_k;
%         Ldl2(i_k) = Ldl2(i_k) + sum(r_k);
%         
%         delta = 0.001;
%         for i_th = 1:8
%             theta_d = theta;
%             theta_d(i_th, i_k) = theta_d(i_th, i_k) + delta;
%             W_kd = update_Wk(ni_hi, nj_hi, ni_lo(i_k), nj_lo(i_k), theta_d(:,i_k), gamma(i_k));
%             r_kd = Y_k - lambda1_k.*(W_kd*x) - lambda2_k;
%             L_kd = sum(r_kd.^2);
%             Ldtheta(i_th, i_k) = (L_kd - L_k) / delta;
%         end
%             
%     end
%     Ldtl = [Ldtheta(:); Ldl1; Ldl2];
% end
% if nargout > 2
%     H = [];
% end
% 
% return
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [L, r] = likelihood_full(x, theta, W, Y, alpha, beta, ni_hi, nj_hi, ni_lo, nj_lo)
% %Construct La and Lb from photometric parameters in theta
% KM = size(Y,1);
% K = length(ni_lo);
% 
% 
% r = Y'-x*W;
% Dx = abs(diff_x(x, ni_hi, nj_hi));
% 
% L = sum(r.^2) + beta*huber_prior(Dx,alpha);
% 
% return

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function e = likelihood_theta_grad(theta, x, W, Y, alpha, beta,ni_hi,nj_hi, ni_lo, nj_lo)
% 
% K = length(ni_lo);
% theta = reshape(theta, [], K);
% e = likelihood_full(x, theta, W, Y, alpha, beta, ni_hi, nj_hi, ni_lo, nj_lo);
% 
% return

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function g = likelihood_x_theta_grad(x_theta, W, Y, alpha, beta, ni_hi, nj_hi, ni_lo, nj_lo)
% 
% K = length(ni_lo);
% num_hi_px = ni_hi*nj_hi;
% x = x_theta(1:num_hi_px);
% %theta = reshape(x_theta(num_hi_px:end), [], K);
% 
% r = Y'-x*W;
% L_dx = W*r';
% L_dx = -2.*L_dx' + elle_eval_huber_grad(x,alpha,beta,ni_hi,nj_hi);
% 
% L_dw = -2
% return
