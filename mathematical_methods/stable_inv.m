function Y = stable_inv(A,B)
% Embedded MATLAB Library function.

% $Revision: 773 $  $Date: 2010-02-02 12:32:07 +0000 (Tue, 02 Feb 2010) $
% Copyright 2002-2005 The MathWorks, Inc.

eml_assert(nargin >= 2, 'Not enough input arguments.');
if isscalar(A)
    Y = B ./ A;
    return
end
eml_assert(isfloat(A), ['Operation ''mldivide'' is not defined for values of class ''' class(A) '''.']);
eml_assert(isfloat(B), ['Operation ''mldivide'' is not defined for values of class ''' class(B) '''.']);
eml_assert(size(B,1) == size(A,1), 'Matrix dimensions must agree.');
eml_assert(ndims(A) == 2 && ndims(B) == 2, 'Input arguments must be 2-D.');
if isempty(A) || isempty(B)
    if isa(A,'single')
        Y = zeros(size(A,2),size(B,2),'single');
    else
        Y = zeros(size(A,2),size(B,2),class(B));
    end
elseif size(A,1) == size(A,2)
    Y = LUSolve(A,B);
else
    Y = QRSolve(A,B);
end

%--------------------------------------------------------------------------

function Y = LUSolve(A,B)
% Solve square A*Y = B via LU decomposition.
if isa(A,'single')
    cls = 'single';
else
    cls = class(B);
end
[X,pivot] = eml_lu(A);
n = size(A,2);
for k = 1:n
    if X(k,k) == 0
        warning('Matrix is singular to working precision.');
        break
    end
end
if (~isreal(A) && isreal(B))
    Y = cast(complex(B(pivot,:)),cls);
else
    Y = cast(B(pivot,:),cls);
end
for k = 1:size(B,2)
    % Solve L*Y = B
    for j = 1:n
        for i = j+1:n
            Y(i,k) = Y(i,k) - Y(j,k)*X(i,j);
        end
    end
    % Solve U*X = Y;
    for j = n:-1:1
        Y(j,k) = Y(j,k) / X(j,j);
        for i = 1:j-1
            Y(i,k) = Y(i,k) - Y(j,k)*X(i,j);
        end
    end
end

%--------------------------------------------------------------------------

function Y = QRSolve(A,B)
% Least squares via QR decomposition.
m = size(A,1);
n = size(A,2);
nb = size(B,2);
mn = min(m,n);
if isa(A,'single')
    cls = 'single';
else
    cls = class(B);
end
[A,tau,jpvt] = eml_zlaqp2(A);
% rankR = rank(R)
tol = max(m,n) * (abs(real(A(1)))+abs(imag(A(1)))) * eps(class(A));
rankR = 0;
for k = 1:mn
    if abs(real(A(k,k))) + abs(imag(A(k,k))) <= tol
        warning(('Rank deficient, rank = %d,  tol = %13.4e.'),rankR,tol);
        break
    end
    rankR = rankR + 1;
end
% Allocate Y
if isreal(A) && isreal(B)
    Y = zeros(n,nb,cls);
else
    Y = complex(zeros(n,nb,cls));
end
% B = Q'*B
if isreal(A) || ~isreal(B)
    for j = 1:mn
        tauj = conj(tau(j));
        if tauj ~= 0
            for k = 1:nb
                wj = cast(B(j,k),cls);
                for i = j+1:m
                    wj = wj + conj(A(i,j))*B(i,k);
                end
                wj = tauj*wj;
                if wj ~= 0
                    B(j,k) = B(j,k) - wj;
                    for i = j+1:m
                        B(i,k) = B(i,k) - A(i,j)*wj;
                    end
                end
            end
        end    
    end    
else
    % Same algorithm but with a complex copy of B.
    CB = complex(cast(B,cls));
    for j = 1:mn
        tauj = conj(tau(j));
        if tauj ~= 0
            for k = 1:nb
                wj = CB(j,k);
                for i = j+1:m
                    wj = wj + conj(A(i,j))*CB(i,k);
                end
                wj = tauj*wj;
                if wj ~= 0
                    CB(j,k) = CB(j,k) - wj;
                    for i = j+1:m
                        CB(i,k) = CB(i,k) - A(i,j)*wj;
                    end
                end
            end
        end    
    end
end
% B(1:rankR,:) = R(1:rankR,1:rankR)\B(1:rankR,:); 
% Y(jpvt,:) = [B; zeros(mb-rankR,nb)];
% We jump around a little bit here to avoid a temporary.
for k = 1:nb
    for i = 1:rankR
        if isreal(B) && ~isreal(A)
            Y(jpvt(i),k) = CB(i,k);
        else
            Y(jpvt(i),k) = cast(B(i,k),cls);
        end
    end
    for j = rankR:-1:1
        pj = jpvt(j);
        Y(pj,k) = Y(pj,k) / A(j,j);
        for i = 1:j-1
            Y(jpvt(i),k) = Y(jpvt(i),k) - Y(pj,k)*A(i,j);
        end
    end
end

%--------------------------------------------------------------------------
