function [consensus sensitivities specificities iter] = staple_consensus(D, varargin)
%STAPLE_CONSENSUS *Insert a one line summary here*
%   [] = staple_consensus(varargin)
%
% STAPLE_CONSENSUS uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments: 
%   D, N*R matrix where each column represents the label of the r-th rater 
%
% Optional Arguments:
%   sensitivities - Initial sensitivities to use for each rater. If given
%   as empty (or not set), default is 0.99999 for each rater. Alternatively,
%   if given a single sensitivity this will be used for all raters, or a
%   1*R vector can be used to set different values for each rater
%
%   specificities - Initial sensitivities to use for each rater (same
%   options available as with sensitivities)
%
%   prior - If set, determine a prior probability of each pixels label being 1
%   (default 0.5 at every pixel)
%
%   max_iterations - Does what it says on the tin...
%
%   tolerance - Minimum change in consensus weights needed to terminate
%   algorithm (assuming max_iterations aren't reached first)
%
% Outputs:
%   consensus - N*1 vector of the final consensus labels at each pixel
%
%   sensitivities - 1*R vector of the final sensitivities of each rater
%
%   specificities - 1*R vector of the final specificities of each rater
%
% Example:
%
% Notes: This function is an implementation of the STAPLE algorithm given
% in the Warfield 2004 paper (IEEE TMI, vol 23, no 7). However this
% algorithm should work for both binary and continuous (on the range 0->1)
% input labels
%
% See also:
%
% Created: 16-Apr-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode % no mandatory arguments
    'sensitivities', [],... %the optional arguments
    'specificities', [],...
    'prior', [], ...
    'tolerance', 1e-4,...
    'max_iterations', 1e3);

[num_pixels num_raters] = size(D);

%Set initial sensitivities based on input arguments
if isempty(args.sensitivities)
    %if empty input, set default sensitivity at 0.99 for all raters
    sensitivities = 0.99999*ones(1, num_raters);
elseif numel(args.sensitivities) == 1
    %equal initial sensitivity supplied for all raters
    sensitivities = args.sensitivities*ones(1, num_raters);
elseif length(args.sensitivities(:)) == num_raters
    %separate initial sensitivity supplied for all raters
    sensitivities = args.sensitivities(:)';
else
    error(['Incorrect number of initial sensitivities given. Must either be empty, 1 or ', num2str(num_raters)]);
end

%Set initial specificities based on input arguments
if isempty(args.specificities)
    %if empty input, set default specifity at 0.99 for all raters
    specificities = 0.99999*ones(num_raters, 1);
elseif numel(args.specificities) == 1
    %equal initial specifity supplied for all raters
    specificities = args.specificities*ones(num_raters, 1);
elseif length(args.specificities(:)) == num_raters
    %separate initial specifity supplied for all raters
    specificities = args.specificities(:)';
else
    error(['Incorrect number of initial specificities given. Must either be empty, 1 or ', num2str(num_raters)]);
end

%set default prior
if isempty(args.prior)
    prior = 0.5*ones(num_pixels, 1);
elseif length(args.prior(:)) == num_pixels
    prior = args.prior(:);
else
    error('If specified the prior must have the same number of elements as pixels in D');
end

%Compute inital label weights as mean of rater labels
W_old = mean(D, 2);

%Set iteration counter and intial change
iter = 1;
change = Inf;

while iter <= args.max_iterations && change > args.tolerance
    
    %E-step: compute label weighting at each pixel
    pos_prob = prior;
    neg_prob = 1 - prior;
    for r = 1:num_raters
        pos_prob = pos_prob .* (sensitivities(r)*D(:,r) + (1-sensitivities(r))*(1-D(:,r)));
        neg_prob = neg_prob .* (specificities(r)*(1-D(:,r)) + (1-specificities(r))*D(:,r));
    end
    
    W_new = pos_prob ./ (pos_prob + neg_prob);
        
    %M-step: compute specifity, sensitivity that maximises the expectation
    %of the data given the new label ratings
    W_sum = sum(W_new);
    for r = 1:num_raters
        sensitivities(r) = sum(D(:,r) .* W_new) / W_sum;
        specificities(r) = sum((1-D(:,r)) .* (1-W_new)) / (num_pixels - W_sum);
    end
    
    %Compute mean absolute change in label weightings
    change = mean(abs(W_old - W_new));
    
    %Save new W into old
    W_old = W_new;
    
    %Increment iteration counter
    iter = iter+1;
end
iter = iter-1;
consensus = W_new;    
    

