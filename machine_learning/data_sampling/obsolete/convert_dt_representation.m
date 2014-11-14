function [dt_samples] = convert_dt_representation(dt_samples, varargin)
%CONVERT_DT_REPRESENTATION *Insert a one line summary here*
%   [] = convert_dt_representation(varargin)
%
% CONVERT_DT_REPRESENTATION uses the U_PACKARGS interface function
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
args = u_packargs(varargin, '0',...
    'feature_type', 'all',...
    'do_max', 0,...
    'win_size', 3,...
    'pca', []);
clear varargin;

if args.do_max
    %get the maximum response across orientations
    dt_samples = max(dt_samples, [], 3);
end

switch args.feature_type
    case 'all'
%         discard_dims = repmat(1:6, 1, size(dt_samples,4)) + 84*kron(0:size(dt_samples,4)-1, ones(1,6));
        dt_samples = [abs(dt_samples(:,:)) angle(dt_samples(:,:))];
%         dt_samples(:, discard_dims) = [];

    case 'mag'
%         discard_dims = repmat(1:6, 1, size(dt_samples,4)) + 84*kron(0:size(dt_samples,4)-1, ones(1,6));
        dt_samples = abs(dt_samples(:,:));
%         dt_samples(:, discard_dims) = [];

    case 'phase'
        dt_samples = angle(dt_samples(:,:));
        
    case 'real'
        dt_samples = real(dt_samples(:,:));

    case 'complex'
        dt_samples = dt_samples(:,:);
        
    case 'conj'
        %Change the sign of samples with negative imaginary part -
        %effectively implemeneting a reflection about the real axis
        fold_idx = imag(dt_samples) < 0;
        dt_samples(fold_idx) = conj(dt_samples(fold_idx));
        dt_samples = [abs(dt_samples(:,:)) angle(dt_samples(:,:))];
        
    case 'ilp'
        %Compute abs of samples
        dt_abs = abs(dt_samples(:,:));

        %Take copy of samples
        ilp_samples = dt_samples;
        
        %Loop through levels
        num_levels = size(dt_samples, 4);
        for lev = 2:num_levels
            %Subtract double the arg of the next coarser level
            ilp_samples(:,:,:,lev) = ...
                dt_samples(:,:,:,lev-1) .* conj(dt_samples(:,:,:,lev).^2);
        end
        
        %Apply conjugate transformation - but reflect in across imag axis
        %fold_idx = real(ilp_samples) < 0;
        %ilp_samples(fold_idx) = complex(-real(ilp_samples(fold_idx)), imag(ilp_samples(fold_idx)));
        %fold_idx = imag(ilp_samples) < 0;
        %ilp_samples(fold_idx) = conj(ilp_samples(fold_idx));
        
        %Combine original DT magnitudes and ILP phases
        dt_samples = [dt_abs angle(ilp_samples(:,:))];
        
    case 'icp'
        %Compute abs of samples
        dt_abs = abs(dt_samples(:,:));

        %Take copy of samples
        icp_samples = dt_samples;
        
        %Workout centre pixel dimension
        cp = (args.win_size^2 + 1) / 2; 
        
        %Loop through the positions in the window
        for pos = 1:args.win_size^2
            if pos == cp; continue; end
            %Subtract the arg of the centre pixel from each surrounding
            %position
            icp_samples(:,pos,:,:) = ...
                dt_samples(:,pos,:,:) .* conj(dt_samples(:,cp,:,:));
        end
        
%         %Apply conjugate transformation
%         fold_idx = imag(icp_samples) < 0;
%         icp_samples(fold_idx) = conj(icp_samples(fold_idx));
        
        %Combine original DT magnitudes and ICP phases
        dt_samples = [dt_abs angle(icp_samples(:,:))];

    otherwise
        error(['Dual tree feature type ', args.feature_type, ' not recognized']);
end 

if ~isempty(args.pca)
    %Transform sample using PCA modes
    dt_samples = bsxfun(@minus, dt_samples, args.pca.mean)*args.pca.modes;
end
