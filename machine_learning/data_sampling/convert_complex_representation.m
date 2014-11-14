function [converted_samples] = convert_complex_representation(complex_samples, feature_type, win_size)
%CONVERT_COMPLEX_REPRESENTATION *Insert a one line summary here*
%   [] = convert_dt_representation(varargin)
%
% CONVERT_COMPLEX_REPRESENTATION uses the U_PACKARGS interface function
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
% Created: 22-Feb-2011
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
if nargin < 2
    win_size = 1;
end

switch feature_type
    case 'all'
        converted_samples = [abs(complex_samples(:,:)) angle(complex_samples(:,:))];
        
    case 'real_imag'
        converted_samples = [real(complex_samples(:,:)) imag(complex_samples(:,:))];
        
    case 'real_abs_imag'
        converted_samples = [real(complex_samples(:,:)) abs(imag(complex_samples(:,:)))];

    case 'mag'
        converted_samples = abs(complex_samples(:,:));

    case 'phase'
        converted_samples = angle(complex_samples(:,:));
        
    case 'real'
        converted_samples = real(complex_samples(:,:));
        
    case 'imag'
        converted_samples = imag(complex_samples(:,:));

    case 'complex'
        converted_samples = complex_samples(:,:);
        
    case 'conj'
        %Change the sign of samples with negative imaginary part -
        %effectively implemeneting a reflection about the real axis
        fold_idx = imag(complex_samples) < 0;
        complex_samples(fold_idx) = conj(complex_samples(fold_idx));
        converted_samples = [abs(complex_samples(:,:)) angle(complex_samples(:,:))];
        
    case 'ilp'
        %Compute abs of samples
        dt_abs = abs(complex_samples(:,:));

        %Take copy of samples
        ilp_samples = complex_samples;
        
        %Loop through levels
        num_levels = size(complex_samples, 4);
        for lev = 2:num_levels
            %Subtract double the arg of the next coarser level
            ilp_samples(:,:,:,lev) = ...
                complex_samples(:,:,:,lev-1) .* conj(complex_samples(:,:,:,lev).^2);
        end
        
        %Apply conjugate transformation - but reflect in across imag axis
        %fold_idx = real(ilp_samples) < 0;
        %ilp_samples(fold_idx) = complex(-real(ilp_samples(fold_idx)), imag(ilp_samples(fold_idx)));
        %fold_idx = imag(ilp_samples) < 0;
        %ilp_samples(fold_idx) = conj(ilp_samples(fold_idx));
        
        %Combine original DT magnitudes and ILP phases
        converted_samples = [dt_abs angle(ilp_samples(:,:))];
        
    case 'icp'
        %Compute abs of samples
        dt_abs = abs(complex_samples(:,:));

        %Take copy of samples
        icp_samples = complex_samples;
        
        %Workout centre pixel dimension
        cp = (win_size^2 + 1) / 2; 
        
        %Loop through the positions in the window
        for pos = 1:win_size^2
            if pos == cp; continue; end
            %Subtract the arg of the centre pixel from each surrounding
            %position
            icp_samples(:,pos,:,:) = ...
                complex_samples(:,pos,:,:) .* conj(complex_samples(:,cp,:,:));
        end
        
        %Apply conjugate transformation
        fold_idx = imag(icp_samples) < 0;
        icp_samples(fold_idx) = conj(icp_samples(fold_idx));
        
        %Combine original DT magnitudes and ICP phases
        converted_samples = [dt_abs angle(icp_samples(:,:))];

    otherwise
        error(['Dual tree feature type ', feature_type, ' not recognized']);
end
