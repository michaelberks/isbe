function [o,gtruth] = synthdata_demo(gtruth, K, transform, gamma, sigma_noise, zm)

%
if ~exist('gtruth', 'var') || isempty(gtruth)
    % westconcordorthophoto.png is a 366x364-pixel greyscale image that comes
    % with Matlab's image processing toolbox. If you don't have it, try
    % substituting any other greyscale image.
    gtruth = imread('westconcordorthophoto.png'); % Load image and scale intensities to [-0.5,0.5]).
end
gtruth = double(gtruth)./255 - 0.5; 

% Create K low-res images using a built-in MATLAB image as ground truth.
if ~exist('K', 'var') || isempty(K)
    K = 5; 
end      % Number of low-res images to generate.

if ~exist('transform', 'var') || isempty(transform)
    transform = 'homography';
end

if ~exist('gamma', 'var') || isempty(gamma)
    gamma = 0.4;     % standard deviation of Gaussian point-spread function for the low-res images.
end

if ~exist('sigma_noise', 'var') || isempty(sigma_noise)
    sigma_noise = 5/255; % standard deviation of noise on the low-res images.
end

if ~exist('zm', 'var') || isempty(zm)
    zm = 2;        % zoom factor for the low-res images.
end

% high-res image size.
[ni_hi, nj_hi] = size(gtruth); 

% This will be my main super-res data structure.
o = struct([]); 
for i = 1:K
    o(i).g = gamma;    % standard deviation of Gaussian point-spread function.
    o(i).n = sigma_noise;  % standard deviation of noise.
    o(i).ni = floor(ni_hi*0.9/zm);     % low-res image size (vertical).(biv)
    o(i).nj = floor(nj_hi*0.9/zm);     % low-res image size (horizontal).(bih)
    o(i).l1 = 1 + randn(1)*0.1;% gsamp(1,0.1^2,1);       % multiplicative photometric parameter (lambda_alpha)
    o(i).l2 = randn(1)*10/255; %gsamp(0,(10/255)^2,1);  % additive photometric parameter (lambda_beta)
    
    switch transform
        case 'homography'
            % Invent a homography using function below.
            o(i).H = local_rand_homography(o(i).ni, o(i).nj, ni_hi, nj_hi, zm); 
        case 'rotation'
            o(i).H = local_rand_rigid(zm);
        case 'translation'
            o(i).H = local_rand_trans(zm);
        case 'identity'
            o(i).H = identity_trans(zm);
    end
    
    % Now generate a low-resolution image using the "gtruth" image and the
    % information we've just put into this structure:
    H = o(i).H';
    theta = H(1:8)';
    lambda = [o(i).l1; o(i).l2];
    
    [o(i).im,o(i).noise,o(i).orig] = ...
        makeLR(gtruth, o(i).ni, o(i).nj, theta, lambda, sigma_noise, gamma);

end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function H = local_rand_homography(ni_lo, nj_lo, ni_hi, nj_hi, zm)
% This creates a random projective homography from a low-res image of size
% smv*smh to a high-res image of size biv*bih, with zoom factor zm.
%
T = diag([1/(nj_lo-1),1/(ni_lo-1),1]); % scaling
shift_lo = [eye(2),-0.5*[nj_lo+1;ni_lo+1];0,0,1]; % shift low-res image coords to be zero-centered.
shift_hi = [eye(2),0.5*[nj_hi+1;ni_hi+1];0,0,1];  % shift high-res image coords to be zero-centered.

P1 = [-0.5,0.5,0.5,-0.5;-0.5,-0.5,0.5,0.5;ones(1,4)];
P2 = [0.15*(rand(2,4)-0.5);zeros(1,4)]+P1;
stack = @(p,q) [zeros(1,3), -p(3)*q', p(2)*q'; p(3)*q', zeros(1,3), -p(1)*q'];
S = []; 
for i = 1:4, 
    S = [S;stack(P1(:,i),P2(:,i))];%#ok 
end 
H = shift_hi*diag([zm,zm,1])*inv(T)*reshape(null(S),3,3)'*T*shift_lo; %#ok
H = H ./ H(9);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function H = local_rand_rigid(zm)
theta = pi*(10*rand - 5) / 180;
xt = zm*(20*rand-10);
yt = zm*(20*rand-10);
H = [...
    zm*cos(theta) zm*sin(theta) xt;
   -zm*sin(theta) zm*cos(theta) yt;
    0            0            1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function H = local_rand_trans(zm)

xt = zm*(20*rand-10);
yt = zm*(20*rand-10);
H = [...
    zm 0 xt;
    0 zm yt;
    0 0  1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function H = identity_trans(zm)
H = [...
    zm 0 0;
    0 zm 0;
    0 0  1];




