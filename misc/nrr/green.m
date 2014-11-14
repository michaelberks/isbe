function [values] = green(xsamp,ysamp,qx0,qy0,params)
%
% Calculates the set of greens functions G( (xsamp,ysamp), (qx0,qy0) )
%
% The type of Greens function is specified in params.green
% Other parameters, like gaussian width, are also in params
%
% If there are nsamp sample points and nknots knots (qx0,qy0)
% then values is size (nsamp x nknots)


types = {'biharmCPS' 'triharmCPS' 'biharmTPS' 'triharmTPS' 'gaussian'};

if(nargin==0)
	values = types;
else
	if( (length(qx0(:))~=length(qy0(:))) | (length(xsamp(:))~=length(ysamp(:))) | isempty(xsamp) )
		values = [];
	else
		% Knot points - 1 x nknots
		qx0 = qx0(:)';
		qy0 = qy0(:)';
		% Sample points - nsamp x 1
		xsamp = xsamp(:);
		ysamp = ysamp(:);

		nknots = length(qx0);
		nsamp = length(xsamp);

		%%%%%%%%%%%%% SWITCH BIT %%%%%%%%%%%%%%%%%%%%%%%%%%
		if(~isfield(params,'green'))
			params.green = char(types(1));
		end

		switch params.green
			case 'biharmTPS'
				% This is NOT the same as the thinplate spline interpolant
				% It is just the Thinplate Greens function!
				for knot = 1:nknots
					p = xsamp.^2 + ysamp.^2;
					q = qx0(knot)^2 + qy0(knot)^2;
					qp = 2*(qx0(knot) * xsamp + qy0(knot) * ysamp);
					pq2 = p - qp + q;
					where = find((pq2==0));
					where2 = find((pq2~=0));					
					values(where,knot) = 0;
					values(where2,knot) = pq2(where2).*log(pq2(where2));
				end
			case 'triharmTPS'
				for knot = 1:nknots
					p = xsamp.^2 + ysamp.^2;
					q = qx0(knot)^2 + qy0(knot)^2;
					qp = 2*(qx0(knot) * xsamp + qy0(knot) * ysamp);
					pq2 = p - qp + q;
					where = find((pq2==0));
					where2 = find((pq2~=0));					
					values(where,knot) = 0;
					values(where2,knot) = (pq2(where2).^2).*log(pq2(where2));
				end
			case 'biharmCPS'
				for  knot = 1:nknots
					p = xsamp.^2 + ysamp.^2;
					q = qx0(knot)^2 + qy0(knot)^2;
					qp = 2*(qx0(knot) * xsamp + qy0(knot) * ysamp);
					PQ2 =  p*q - qp + 1;
					pq2 = p - qp + q;
					where = find((pq2==0));
					where2 = find((pq2~=0)&(PQ2~=0));
					where3 = find((pq2~=0)&(PQ2==0));					
					values(where2,knot) = PQ2(where2) - pq2(where2) - pq2(where2).*log( PQ2(where2)./pq2(where2) );
					values(where,knot) = PQ2(where);
					values(where3,knot) = 0;
					% Points outside unit circle
					where4 = find(p>1);
					values(where4,knot) = 0;
				end
			case 'gaussian'
				% Default value of sigma
				if(~isfield(params,'sigma'))
					params.sigma = 0.1;
				end
				for knot = 1:nknots
					values(:,knot) = exp( - ( (xsamp-qx0(knot)).^2 + (ysamp-qy0(knot)).^2 )/(2*params.sigma*params.sigma));
				end
			case 'triharmCPS'
				for  knot = 1:nknots
					p = xsamp.^2 + ysamp.^2;
					q = qx0(knot).^2 + qy0(knot).^2;
					qp = 2*(qx0(knot) * xsamp + qy0(knot) * ysamp);
					PQ =  sqrt(p*q - qp + 1);
					pq = sqrt(p - qp + q);
					where = find(pq==0);
					where2 = find((pq~=0)&(PQ~=0));
					where3 = find((pq~=0)&(PQ==0));					
					values(where2,knot) = 0.25*(PQ(where2).^4) - (pq(where2).^2).*(PQ(where2).^2) + 0.75*(pq(where2).^4) + (pq(where2).^4).*log(PQ(where2)./pq(where2));
					values(where,knot) = 0.25*(PQ(where).^4);
					values(where3,knot) = 0;
					% Points outside unit circle
					where4 = find(p>1);
					values(where4,knot) = 0;
				end
		%%%%%%%%%% END OF SWITCH %%%%%%%%%%%%%%%%%%%%%%%
		end
	end
end


