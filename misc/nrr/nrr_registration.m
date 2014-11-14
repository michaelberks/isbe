function [data,params] = nrr_registration(data,params)
% params.nrr_type : integer, different optimisation schemes
% eg 1 : all knots, one at a time, random order
% 
% Other parameters = type specific
% eg. params.nloops, params.weights, blah blah
%
% Also params.nrr_options for fminsearch - may be more than one lot
%
%
%
% Optimiser: params.nrrOJ : objective function type
%

if(~isfield(data,'q0'))
	data.q0 = [];
	data.q1 = [];	
elseif(isempty(data.q0))
	data.q0 = [];
	data.q1 = [];	
else

if(~isfield(params,'showconv'))
	params.showconv = 0;
end

% Keeps running count of OJ value as function of number of function evaluations
if(params.showconv==1)
	fcount = [];
	OJvalue = [];
end



%%%%%%%%%%%%%%%%%%% MAIN TYPE LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	switch params.nrr_type
		case 1
		% Run over all knotpoints, random order,nloops times, one knot at a time
		if(~isfield(data,'q1'))
			data.q1 = data.q0;
		elseif(isempty(data.q1))
			data.q1 = data.q0;
		end
		nknots = size(data.q0,2);
		for i=1:params.nloops_nrr
			% knots in random order - equally weighted
			pset = randperm(nknots)
			for j=pset
				whichvar = j;
				whichfix = [1:nknots].*([1:nknots]~=j);
				where = find(whichfix~=0);
				whichfix = whichfix(where);
				[data.q1(:,whichvar(:)),OJbest,dummy,output] = fminsearch(@eval_nrr_objective_function,data.q1(:,whichvar(:)),params.nrr_options,data.q1(:,whichfix(:)),whichvar,whichfix,data.q0(1,:),data.q0(2,:),data,params);
				if(params.showconv==1)
					if(isempty(fcount))
						fcount = output.funcCount;
						OJvalue = OJbest;
					else
						OJvalue = [OJvalue OJbest];	
						fcount(length(fcount)+1) = fcount(length(fcount)) + output.funcCount;
					end
				end				
			end
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		case 2
		% Weighted random order over all knotpoints, one at a time
		if(~isfield(data,'q1'))
			data.q1 = data.q0;
		elseif(isempty(data.q1))
			data.q1 = data.q0;
		end
		nknots = size(data.q0,2);
		weights = get_weights(data,params);
		% Random permutation, weighted by weights!
		nset = [];
		for i=1:nknots
			n = round(params.nloops*nknots*weights(i));
			nset(length(nset)+1:length(nset)+n) = i;
		end
		nset = nset(randperm(length(nset)));

		% For animation, draw figure
		if(strcmp(params.anim,'yes'))
			[animation,params] = ani_warp(data,params);
		end

		for j=nset
				whichvar = j;
				whichfix = [1:nknots].*([1:nknots]~=j);
				where = find(whichfix~=0);
				whichfix = whichfix(where);
				[data.q1(:,whichvar(:)),OJbest,dummy,output] = fminsearch(@eval_nrr_objective_function,data.q1(:,whichvar(:)),params.nrr_options,data.q1(:,whichfix(:)),whichvar,whichfix,data.q0(1,:),data.q0(2,:),data,params);
				if(params.showconv==1)
					if(isempty(fcount))
						fcount = output.funcCount;
						OJvalue = OJbest;
					else
						OJvalue = [OJvalue OJbest];	
						fcount(length(fcount)+1) = fcount(length(fcount)) + output.funcCount;
					end
				end

				% Animation
				if(strcmp(params.anim,'yes'))
					[animation,params] = ani_warp(data,params,animation);
				end

				
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		case 3
		% One by one newer points, perturb older

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		case 4
		% Some parallel stuff

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	otherwise
		data.q1 = [];
	%%%%%% END OF SWITCH %%%%%%%%%%%%%%%	
	end
%%%% END OF MAIN LOOP %%%%%%%%%%%%%%%%%%
end

% data.age??

% Regularisation term in eval_nrr_objective function = parameters in params

% Create the nrr warped reference grid
[xgrid_nrr,ygrid_nrr,dummy,data] = nrr_trans(data.xgrid_affine,data.ygrid_affine,data.q0(1,:),data.q0(2,:),data.q1(1,:),data.q1(2,:),params,data);
data.xgrid_nrr = xgrid_nrr;
data.ygrid_nrr = ygrid_nrr;

% Plot the convergence
if(params.showconv==1)
	figure, plot(fcount,OJvalue,'ko-');
end				

