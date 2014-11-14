function [boostreg] = pt_boosted_reg_train(X, t, varargin)
%PT_BOOSTED_REG_TRAIN Train a boosted regressor from training pairs (X,t)
%   [boostreg] = pt_boosted_reg_train(X, t, varargin)
%
% Inputs:
%      X:      N x D matrix of input data, where each row is a
%              datapoint consisting of d input variables
%
%      t:      N x d vector of class labels for each data point
%
%
% Outputs:
%      boostreg - structure containing regressor parameters
%
%
% Example:
%
% Notes:
%
% See also:
%   TREE_REG_TRAIN, PT_LIN_REG_TRAIN
%
% Created: 10-Mar-2011
% Author: Phil Tresadern 
% Email : philip.tresadern@manchester.ac.uk 
% Phone : +44 (0)161 275 5114 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin,... % the user's input
			 '0', ... % non-strict mode
			 'boost_n_levels', 100, ...
             'boost_weak_learner', 'piecewise_constant', ...
			 'boost_output_type', 'raw', ...
			 'boost_n_bins', 24, ...
			 'boost_shrinkage', 0.05);

% Lose the boost_ prefix for clarity
args = remove_prefix('boost_',args);
         
% check dimensions are consistent
[Nx,Dx] = size(X);
[Nt,Dt] = size(t);
if (Nx~=Nt)
	error('Different number of inputs to outputs');
end
N = Nx;

% % find and move any NaNs in outputs
% % if matrix is very big, deleting rows causes menory upsets
% nan_inds = find(any(isnan(t),2));
% if ~isempty(nan_inds)
% 	t(nan_inds,:) = t(end-length(nan_inds)+1:end,:);
% 	X(nan_inds,:) = X(end-length(nan_inds)+1:end,:);
% 	N = N-length(nan_inds);
% end

nan_inds = find(any(isnan(t),2));
if ~isempty(nan_inds)
	t(nan_inds,:) = [];
	X(nan_inds,:) = [];
	N = size(X,1);
end

% get shrinkage parameter
if ~isscalar(args.shrinkage)
	error('Shrinkage argument must be a scalar');
end
shrinkage = args.shrinkage;

% preprocess the inputs
switch args.weak_learner
	case 'piecewise_constant',
		% compute bin edges such that samples are uniformly distributed
		X_bins = X;
		bins = ceil(args.n_bins*(1:N)/N);
		edge_inds = find(diff(bins));
		bin_edges = zeros(args.n_bins-1,Dx);
		for d = 1:Dx
			[X_sorted,inds] = sort(X(:,d));
			X_bins(inds,d) = bins;
			bin_edges(:,d) = X_sorted(edge_inds);
		end
end

% compute initial residual (i.e. t)
resid = t;

% begin boosting
for L = 1:args.n_levels
	
	% initial error
	err_min = inf;
	
	% search over inputs for best predicting one
	for inp = 1:Dx

		% initialize residual for this input variable
		resid_i = resid;
		
		switch args.weak_learner
			case 'piecewise_constant',
				% initialize output
				y = zeros(1,args.n_bins);
				
				for bin = 1:args.n_bins
					% find which samples are in this bin
					inds = (X_bins(:,inp)==bin);
					if isempty(inds), continue; end
					
					% set prediction to mean value for this bin
					y(bin) = mean(resid(inds));

					% compute corresponding residual
					switch args.output_type
						case 'raw',
							resid_i(inds) = resid_i(inds) - y(bin)*shrinkage;
						case 'angle',
							resid_i(inds) = resid_i(inds) * ...
											exp(complex(0,-angle(y(bin))*shrinkage));
					end
				end
				
				% compute error for this input variable
				switch args.output_type
					case 'raw',
						err = mean(abs(resid_i));
					case 'angle',
						err = 1-abs(mean(resid_i));
				end
				
				% if error lower than others then store parameters
				if (err < err_min)
					best_input = inp;
					best_y = y;
					err_min = err;
				end
				
			otherwise,
				error(['Unknown weak learner type: ',args.weak_learner]);
		end
	end

	% deal with missing values, shrink y's contribution and update residual
	switch args.output_type
		case 'raw',
			best_y(isnan(best_y)) = 0;
			best_y = shrinkage*best_y;
			for bin = 1:args.n_bins
				inds = (X_bins(:,best_input)==bin);
				resid(inds) = resid(inds) - best_y(bin);
			end
		case 'angle',
			best_y(isnan(best_y)) = complex(1,0);
			best_y = exp(complex(0,angle(best_y)*shrinkage));
			for bin = 1:args.n_bins
				inds = (X_bins(:,best_input)==bin);
				resid(inds) = resid(inds) * conj(best_y(bin));
			end
	end
	
	% store best parameters
	levels(L) = struct('input',best_input,'y',best_y);
end

% store arguments along with regressor parameters
boostreg = args;

% bin parameters for each dimension
boostreg.bin_edges = bin_edges;

% parameters for each learner in ensemble
boostreg.levels = levels;
