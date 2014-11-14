function [mu_hat, kappa_hat] = von_mises_mle(sample)
% Compute the maximum liklihood estimates for mu and kappa assuming sample
% comes from Von Mises distribution

% Compute C_ and S_
C_ = mean(cos(sample));
S_ = mean(sin(sample));

%Compute R_ and theta_
R_ = sqrt(C_^2 + S_^2);
theta_ = atan2(S_,C_);

% Now compute the MLE estimates
mu_hat = theta_;

% kappa_hat = A^-1(R_);
kappa_hat = Ainv(R_);




