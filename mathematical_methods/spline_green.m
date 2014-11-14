function [values] = spline_green(src_xy, int_xy, green)
%
% Calculates the set of greens functions G( X, Y )
%
% The type of Greens function is specified in params.green
% Other parameters, like gaussian width, are also in params
%
% If there are nsamp sample points and nknots knots (src_x,src_y)
% then values is size (nsamp x nknots)
if nargin < 3
    green = 'biharmTPS';
end

if isempty(int_xy)
    int_xy = src_xy;
end

%%%%%%%%%%%%% SWITCH BIT %%%%%%%%%%%%%%%%%%%%%%%%%%

switch green
    case 'biharmTPS'
        % This is NOT the same as the thinplate spline interpolant
        % It is just the Thinplate Greens function!

        nc = size(src_xy, 2);
        nx = size(int_xy, 2);

        temp = repmat(1:nx,nc,1);
        points = int_xy(:,temp) - reshape(repmat(src_xy(:),1,nx),2,nx*nc);
        ap2 = sum(points.*points,1); ap2(ap2==0) = 1;
        values= reshape(ap2.*log(ap2),nc,nx)';

   case 'biharmCPS' 
        nc = size(src_xy, 2);
        nx = size(int_xy, 2);

        int_x= int_xy(1,:)';
        int_y= int_xy(2,:)';

        src_x= src_xy(1,:);
        src_y= src_xy(2,:);
        
        p = repmat(int_x.^2 + int_y.^2,1,nc);
        q = repmat(src_x.^2 + src_y.^2,nx,1);
        
        qp= 2*(src_x'*int_x' + src_y'*int_y')';
        pq2= abs(p - qp + q);
        PQ2= p(:,1)*q(1,:) - qp + 1;
        ind= pq2&PQ2;
        values= zeros(size(p));
        values(ind)= PQ2(ind) - pq2(ind) - pq2(ind).*log(PQ2(ind)./pq2(ind));
        values(~pq2)= PQ2(~pq2);
        values(pq2&~PQ2)= 0;
        values(p>1)= 0;
        %%%%%%%%%% END OF SWITCH %%%%%%%%%%%%%%%%%%%%%%%
end

values = values';