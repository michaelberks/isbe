function [bool] = mb_issymmetric(M)
%
% mb_issymmetric(M)
%
% Tests if M is symmetric

if size(M,1)~=size(M,2),
  bool = false; return
else
  tol = eps^(2/3);
end

bool = max(max(abs(M-M'))) <= tol*max(abs(M(:)));
