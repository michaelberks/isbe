function [c] = mb_sample_based_on_probs(varargin)
%
% MB_SAMPLE_BASED_ON_PROBS Sample based on probability or frequency count.
%
% The set of events I = [I1, I2, ..., IN], has associated probabilities
% P = [P1, P2, ..., PN]. This function returns the index c to one of the events
% in I according to the pdf defined by P. If P is actually a set of frequency counts
% P = [F1, F2, ..., FN], then these are normalised and treated as probabilities.
% There is a flag that can be set to ensure that erroneous probability vectors
% are not misintepreted as frequency counts.
%
% MB_SAMPLE_BASED_ON_PROBS uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory:
%
% 'P'
%		- a vector of either probabilities (where sum(P)==1) or frequency counts
%		(where sum(P)>0 & sum(P)~=1).
%
% 'Strict'
%		- a string which indicates whether P must be probabilities (i.e. sum(P) must
%		sum to unity), or whether it can be considered as a set of frequency counts.
%		Must be either 'strict' or 'notstrict'. If 'strict', then the P vector must sum
%		to unity (i.e. be able to be treated as probabilities directly), else an error is
%		reported. If 'notstrict', then the P vector does not have to sum to unity, and
%		if it does not, then the values in P are assumed to be frequency counts and
%		are normalised to unity, and then treated as probabilities.
%
% MB_SAMPLE_BASED_ON_PROBS returns C, which is the index to the event selected.

% Method is the one suggested by Knuth, in Vol 2 of the Art of Computer Programming,
% 3rd edition, in section 3.4.1. 

% pack the args
args = u_packargs(varargin, ... % what the user sent us
			'notstrict',...
			{'Strict', 'P'},...
			'Strict', [],...
			'P', []);


% check that the strict argument is specified properly
if ~(strcmp(args.Strict, 'strict') || strcmp(args.Strict, 'notstrict'))
	error('The Strict argument is not specified correctly; it must be ''strict'' or ''notstrict''')
end

% check that probs is not zero
if sum(args.P) == 0
	error('All the probs were zero, so we can''t sample')
end

% check whether to treat probs as being frequency counts
if ~fp_compare(sum(args.P), 1)
	if strcmp(args.Strict, 'notstrict')% treat probs as being a frequency count; normalise to get probabilities
		args.P = args.P ./ sum(args.P);
	else
		disp('ERROR: The probabilities do not sum to unity, check to make sure your probability vector was computed correctly,')
		error('or specify the ''notstrict'' option to treat the vector as frequency counts')
	end
end

% select a value from a uniform dist in range 0 to 1
u = rand(1);
if 0 <= u && u < args.P(1)
	c = 1;
	return;
end
for i = 2 : length(args.P)
	if args.P(i-1) <= u && u < sum(args.P(1:i))
		c = i;
		break;
	end
end

if isempty(c)
	error('There is an algorithmic problem in MB_SAMPLE_BASED_ON_PROBS')
end
