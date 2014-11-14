% PINV_COV	Compute the pseudo-inverse of a covariance matrix
%
% X = pinv_cov(S)
%
% This function computes the pseudo-inverse of a covariance matrix
% (i.e. a real, symmetric matrix with non-negative diagonals).
% The pseudo-inverse is generally considered to be a robust way of
% computing the inverse of covariance matrices, which often appear
% to be ill-conditioned. However, the method is not foolproof, and
% you may have to take additional steps (e.g. adding variance to the
% diagonal).
%
% This function is implemented using a MEX function (pinv_cov.c), and
% is generally faster that Matlab's PINV function. Experience suggests
% that, in a real-world situation, using PINV_COV rather than PINV can
% double the speed of your code (assuming the main thing your code does is
% compute inverses, e.g. conditioning a GMM).
%
% Because speed was the goal rather than idiot-proof code, very little error
% checking is performed. DO NOT send matrices that are not covariance
% matrices to this function; doing so will probably result in a
% segmentation error and you will have to restart Matlab. Or, even worse,
% you will get numerically incorrect answers!
%
% To obtain the best performance from this function, ensure that your matrices
% are positive definite (i.e. do PCA!). See the C source file for more details.
