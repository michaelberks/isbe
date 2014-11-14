function [xout,yout,E,data] = nrr_trans(xin,yin,qx0,qy0,qx1,qy1,params,data)
% Does a nrr transformation on the points xin yin
%
% (qx0,qy0) are the initial knotpoint positions, (qx1,qy1) the final
% params.green specifies the type of Greens function to use

nknots = length(qx0(:));
% make all knotpoint positions size (nknots x 1)
qx0 = qx0(:);
qy0 = qy0(:);
qx1 = qx1(:);
qy1 = qy1(:);
if( (size(qx0)~=size(qy0)) | (size(qx0)~=size(qx1)) |(size(qx0)~=size(qy1)) )
	xout = [];
	yout = [];
else
	if(isempty(qx0))
		xout = xin;
		yout = yin;
	else	
        % Only need to calculate this once!
        if(~isfield(data,'G'))
		    % Calculate G between all knot point pairs
		    data.G = green(qx0,qy0,qx0,qy0,params);
			data.Ginv = pinv(data.G);
               
        elseif ((size(qx0,1)) ~= (size(data.G,1)))
            % Calculate G between all knot point pairs
	        data.G = green(qx0,qy0,qx0,qy0,params);
            data.Ginv = pinv(data.G);
        end
    alphax = data.Ginv*(qx1-qx0);
	alphay = data.Ginv*(qy1-qy0);

	% Calculate E
	E = alphax'*data.G*alphax + alphay'*data.G*alphay;	

	% Sampled points
	orig_size = size(xin);
	xin = xin(:);
	yin = yin(:);
	xout = xin;
	yout = yin;
	if(~isempty(xin))
		S = green(xin,yin,qx0,qy0,params);
		xout = xout + S*alphax;
		yout = yout + S*alphay;
	end
	xout = reshape(xout,orig_size);
	yout = reshape(yout,orig_size);
    end
end