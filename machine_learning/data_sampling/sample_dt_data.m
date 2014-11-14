function [dt_samples] = sample_dt_data(dt, rows, cols, varargin)
%SAMPLE_DT_DATA *Insert a one line summary here*
%   [] = sample_dt_data(varargin)
%
% SAMPLE_DT_DATA uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%
% Optional Arguments:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 25-Aug-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin, '0',...
    'feature_shape', 'rect',...
    'levels', 1:5,...
    'do_max', 0,...
    'rotate', 0,...
    'win_size', 3,...
    'interpmethod', 'cubic',...
    'band_fequencies', 1i*[-3 -1; -sqrt(5) -sqrt(5); -1 -3; 1 -3; sqrt(5) -sqrt(5); 3 -1]*pi/2.15,...
    'unwrap_phase', 1,...
    'interp_mag_phase', 0,...
    'correct_phase', 1);

% Now sample the training data
switch lower(args.feature_shape)
    case {'clock'}
        [dt_samples] = sample_dt_data_clock(dt, rows, cols, ...
                            args.win_size, args.rotate, args.levels);

    case {'polar7'}
        [dt_samples] = sample_dt_data_polar7(dt, rows, cols, ...
                            args.rotate);
        
    case {'polar13'}
        [dt_samples] = sample_dt_data_polar13(dt, rows, cols, ...
                            args.rotate);

    case {'rect'}
        [dt_samples] = sample_dt_data_rectangle(dt, rows, cols, ...
                            args);    

    otherwise
        error(['Dual tree feature shape ', args.feature_shape, ' not recognized']);
end

% Done with these now
clear dt rows cols;
