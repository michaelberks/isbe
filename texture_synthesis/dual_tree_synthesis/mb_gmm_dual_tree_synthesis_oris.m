function [synthesised_image, dual_tree, cluster_image] = mb_gmm_dual_tree_synthesis_oris(varargin)
%
%MB_GMM_DUAL_TREE_SYNTHESIS Synthesise an image using an ergodic GMM.
%
% 
%
% MB_GMM_DUAL_TREE_SYNTHESIS uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%
% 'FilledImage'
%   - a binary image that indicates where
%
% Optional Arguments:
%
%
% 'SynthesisMode'
%   - specifies whether synthesis is
%
% Return Value:
%
%   MB_GMM_DUAL_TREE_SYNTHESIS returns SYNTHESISED_IMAGE
%
% Example Usage:
%
%

%
% References:
%

% Unpack the arguments:
args = u_packargs(varargin, '0', ... % strict mode
		  {...  % The mandatory arguments
          'FilledImage',...
          'ModelDir',...
          'ModelName',...
          },... % The optional arguments
          'TargetImage', [],...
          'TargetDualTree', [],...
          'ClusterImage', [], ...
          'SynthesisMode', 'simple', ...
		  'SaveFile', [],...
          'NumLevels', 5,...
          'CutOffLevel', 4,...
          'ForceCluster', 0, ...
          'Plot', 0,...
          'MovieFilename', []);

clear varargin;

display(args);

% Use a proxy for the 'SynthesisMode' parameter, to avoid doing
% expensive string comparisons during run-time.
switch args.SynthesisMode
 case 'simple'
  
 case 'advanced'
  
 otherwise
  error(['The ' args.SynthesisMode ' mode is not supported.']);
end

%check model directory is filesep terminated
if ~strcmp(args.ModelDir(end), filesep)
	args.ModelDir = [args.ModelDir filesep];
end

if ~isempty(args.TargetImage)
    % convert TargetImage to doubles, maintain the pixel value ranges (0-255)
    args.TargetImage = double(args.TargetImage); % TargetImage is now a double
    
    % Calculate the dual_tree for the image to be synthesised
    [dual_tree] = dtwavexfm2(args.TargetImage, args.NumLevels,'near_sym_b','qshift_b');

elseif ~isempty(args.TargetDualTree)
    dual_tree = args.TargetDualTree;
    args = rmfield(args, 'TargetDualTree');
else
    error('Either a target image or target dual_tree must be supplied');
end

num_orientations = size(dual_tree{1}, 3);

% get row/col subscripts of unfilled pixels
[p_rows p_cols] = find(~args.FilledImage);

tic; % Start a record of when we start actually working on the synthesis

%pre-allocate cluster_image for output
cluster_image = cell(args.NumLevels, num_orientations);

for level = args.CutOffLevel:-1:1
    
    
    % Calculate indices for that level and make unique
    new_rows = ceil(p_rows/2^(level-1));
    new_cols = ceil(p_cols/2^(level-1));
    new_idx = sub2ind(size(dual_tree{level}(:,:,1)), new_rows, new_cols);

    %make filled image from new_idx
    filled_image = ones(size(dual_tree{level}(:,:,1)));
    filled_image(new_idx) = 0;
    
    syn_args.PixelOrder = mb_compute_pixel_order(filled_image);
    
    %set force clusters flag
    syn_args.ForceCluster = args.ForceCluster;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Synthesise real bands {1,2,3} 
    if ~isempty(args.MovieFilename)
        syn_args.MovieFilename = ...
            [args.MovieFilename, '_',...
        num2str(level), '_1_1.gif'];
    end
    
    syn_args.PathToTextureGMM = ...
        [args.ModelDir, args.ModelName, '_',...
        num2str(level), '_1_1.mat'];
        
    if exist(syn_args.PathToTextureGMM, 'file')

        display(['Synthesising from ', syn_args.PathToTextureGMM]);
        syn_args.UpperBand1 = real(dual_tree{level}(:,:,1));
        syn_args.UpperBand2 = real(dual_tree{level}(:,:,2));
        syn_args.UpperBand3 = real(dual_tree{level}(:,:,3));
        syn_args.LowerBand1 = real(dual_tree{level+1}(:,:,1));
        syn_args.LowerBand2 = real(dual_tree{level+1}(:,:,2));
        syn_args.LowerBand3 = real(dual_tree{level+1}(:,:,3));

        if ~isempty(args.ClusterImage)
            syn_args.ClusterImage = args.ClusterImage{level, 1};
        end

        %Now synthesise the new texture for the band and save back to
        %dual_tree
        [real_123 cluster_image{level, 1}] = ...
                mb_gmm_tex_synthesis_oris(syn_args);
    else
        %don't try and synthesise from a model that doesn't exist!
        warning(['Model ', syn_args.PathToTextureGMM, ' does not exist']); %#ok
        real_123 = {real(dual_tree{level}(:,:,1)),...
                    real(dual_tree{level}(:,:,2)),...
                    real(dual_tree{level}(:,:,3))};
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Synthesise real bands {4,5,6} 
    if ~isempty(args.MovieFilename)
        syn_args.MovieFilename = ...
            [args.MovieFilename, '_',...
        num2str(level), '_4_1.gif'];
    end
    
    syn_args.PathToTextureGMM = ...
        [args.ModelDir, args.ModelName, '_',...
        num2str(level), '_4_1.mat'];
        
    if exist(syn_args.PathToTextureGMM, 'file')

        display(['Synthesising from ', syn_args.PathToTextureGMM]);
        syn_args.UpperBand1 = real(dual_tree{level}(:,:,4));
        syn_args.UpperBand2 = real(dual_tree{level}(:,:,5));
        syn_args.UpperBand3 = real(dual_tree{level}(:,:,6));
        syn_args.LowerBand1 = real(dual_tree{level+1}(:,:,4));
        syn_args.LowerBand2 = real(dual_tree{level+1}(:,:,5));
        syn_args.LowerBand3 = real(dual_tree{level+1}(:,:,6));

        if ~isempty(args.ClusterImage)
            syn_args.ClusterImage = args.ClusterImage{level, 2};
        end

        %Now synthesise the new texture for the band and save back to
        %dual_tree
        [real_456 cluster_image{level, 2}] = ...
                mb_gmm_tex_synthesis_oris(syn_args);
    else
        %don't try and synthesise from a model that doesn't exist!
        warning(['Model ', syn_args.PathToTextureGMM, ' does not exist']); %#ok
        real_456 = {real(dual_tree{level}(:,:,4)),...
                    real(dual_tree{level}(:,:,5)),...
                    real(dual_tree{level}(:,:,6))};
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Synthesise imaginary bands {1,2,3} 
    if ~isempty(args.MovieFilename)
        syn_args.MovieFilename = ...
            [args.MovieFilename, '_',...
        num2str(level), '_1_0.gif'];
    end
    
    syn_args.PathToTextureGMM = ...
        [args.ModelDir, args.ModelName, '_',...
        num2str(level), '_1_0.mat'];
        
    if exist(syn_args.PathToTextureGMM, 'file')

        display(['Synthesising from ', syn_args.PathToTextureGMM]);
        syn_args.UpperBand1 = imag(dual_tree{level}(:,:,1));
        syn_args.UpperBand2 = imag(dual_tree{level}(:,:,2));
        syn_args.UpperBand3 = imag(dual_tree{level}(:,:,3));
        syn_args.LowerBand1 = imag(dual_tree{level+1}(:,:,1));
        syn_args.LowerBand2 = imag(dual_tree{level+1}(:,:,2));
        syn_args.LowerBand3 = imag(dual_tree{level+1}(:,:,3));

        if ~isempty(args.ClusterImage)
            syn_args.ClusterImage = args.ClusterImage{level, 3};
        end

        %Now synthesise the new texture for the band and save back to
        %dual_tree
        [imag_123 cluster_image{level, 3}] = ...
                mb_gmm_tex_synthesis_oris(syn_args);
    else
        %don't try and synthesise from a model that doesn't exist!
        warning(['Model ', syn_args.PathToTextureGMM, ' does not exist']); %#ok
        imag_123 = {imag(dual_tree{level}(:,:,1)),...
                    imag(dual_tree{level}(:,:,2)),...
                    imag(dual_tree{level}(:,:,3))};
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Synthesise imag bands {4,5,6} 
    if ~isempty(args.MovieFilename)
        syn_args.MovieFilename = ...
            [args.MovieFilename, '_',...
        num2str(level), '_4_0.gif'];
    end
    
    syn_args.PathToTextureGMM = ...
        [args.ModelDir, args.ModelName, '_',...
        num2str(level), '_4_0.mat'];
        
    if exist(syn_args.PathToTextureGMM, 'file')

        display(['Synthesising from ', syn_args.PathToTextureGMM]);
        syn_args.UpperBand1 = imag(dual_tree{level}(:,:,4));
        syn_args.UpperBand2 = imag(dual_tree{level}(:,:,5));
        syn_args.UpperBand3 = imag(dual_tree{level}(:,:,6));
        syn_args.LowerBand1 = imag(dual_tree{level+1}(:,:,4));
        syn_args.LowerBand2 = imag(dual_tree{level+1}(:,:,5));
        syn_args.LowerBand3 = imag(dual_tree{level+1}(:,:,6));

        if ~isempty(args.ClusterImage)
            syn_args.ClusterImage = args.ClusterImage{level, 4};
        end

        %Now synthesise the new texture for the band and save back to
        %dual_tree
        [imag_456 cluster_image{level, 4}] = ...
                mb_gmm_tex_synthesis_oris(syn_args);
    else
        %don't try and synthesise from a model that doesn't exist!
        warning(['Model ', syn_args.PathToTextureGMM, ' does not exist']); %#ok
        imag_456 = {imag(dual_tree{level}(:,:,4)),...
                    imag(dual_tree{level}(:,:,5)),...
                    imag(dual_tree{level}(:,:,6))};
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    dual_tree{level}(:,:,1) = complex(real_123{1}, imag_123{1});
    dual_tree{level}(:,:,2) = complex(real_123{2}, imag_123{2});
    dual_tree{level}(:,:,3) = complex(real_123{3}, imag_123{3});
    dual_tree{level}(:,:,4) = complex(real_456{1}, imag_456{1});
    dual_tree{level}(:,:,5) = complex(real_456{2}, imag_456{2});
    dual_tree{level}(:,:,6) = complex(real_456{3}, imag_456{3});
    
    %display new texture if required
    if args.Plot
        for ori = 1:6
            figure; 
            subplot(1,2,1); imagesc(real(dual_tree{level}(:,:,ori)));...
                axis image; colormap(jet(256));
            subplot(1,2,2); imagesc(imag(dual_tree{level}(:,:,ori)));...
                axis image; colormap(jet(256));
            %subplot(2,2,3:4); imagesc(cluster_image{level, ori});...
            %    axis image; colormap(jet(256));
        end
    end
end

% reconstruct image from dual_tree
synthesised_image = uint8(dtwaveifm2(dual_tree,'near_sym_b','qshift_b'));

%get filled image to save with results
filled_image = args.FilledImage; %#ok

% save final image
if ~isempty(args.SaveFile)
    save(args.SaveFile, 'synthesised_image', 'dual_tree', 'cluster_image', 'filled_image');
end

% Report the generation of the synthetic image
disp('Progress: 100 percent complete.');
disp('A modified mamogram region has been synthesised!');

%Main function end