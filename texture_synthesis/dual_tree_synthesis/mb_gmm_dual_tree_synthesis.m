function [synthesised_image, dual_tree, cluster_image] = mb_gmm_dual_tree_synthesis(varargin)
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
          'SynthesisMode', 'patch-wise', ...
		  'SaveFile', [],...
          'NumLevels', 5,...
          'CutOffLevel', 4,...
          'ForceCluster', 0, ...
          'Plot', 0,...
          'ConditionLevels', 1,...
          'MovieFilename', []);

clear varargin;

display(args);

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
    syn_args.SynthesisMode = args.SynthesisMode;
    
    %make sure uint8 mode is turned off since we're synthesising doubles
    %and complex mode is turned on
    syn_args.Uint8Mode = 0;
    
    %set force clusters flag
    syn_args.ForceCluster = args.ForceCluster;
    
    for ori = 1:num_orientations
        
        %Set arguments specific to current band
        
        if ~isempty(args.MovieFilename)
            syn_args.MovieFilename = ...
                [args.MovieFilename, '_',...
            num2str(level), '_', num2str(ori), '.gif'];
        end
        
        syn_args.PathToTextureGMM = ...
            [args.ModelDir, args.ModelName, '_',...
            num2str(level), '_', num2str(ori), '.mat'];
        
        if exist(syn_args.PathToTextureGMM, 'file')
            
            display(['Synthesising from ', syn_args.PathToTextureGMM]);
            syn_args.TargetImage = dual_tree{level}(:,:,ori);
            
            if ~isempty(args.ClusterImage)
                syn_args.ClusterImage = args.ClusterImage{level, ori};
            end
            
            %Now synthesise the new texture for the band and save back to
            %dual_tree
            %Work out whether we're conditioning on the lower level
            if args.ConditionLevels
                %Use mb_gmm_tex_synthesis2
                %Need to give lower dual_tree level as sample image
                if level + 1 == size(dual_tree, 1)
                    syn_args.SampleImage2 = dual_tree{level+1};
                else
                    syn_args.SampleImage2 = dual_tree{level+1}(:,:,ori);
                end
                
                [dual_tree{level}(:,:,ori) cluster_image{level, ori}] = ...
                        mb_gmm_tex_synthesis2_complex(syn_args);
            else
                %Use mb_gmm_tex_synthesis
                [dual_tree{level}(:,:,ori) cluster_image{level, ori}] = ...
                    mb_gmm_tex_synthesis_complex(syn_args);
            end

            %display new texture if required
            if args.Plot
                figure; 
                subplot(2,2,1); imagesc(real(dual_tree{level}(:,:,ori)));...
                    axis image; colormap(jet(256));
                subplot(2,2,2); imagesc(imag(dual_tree{level}(:,:,ori)));...
                    axis image; colormap(jet(256));
                subplot(2,2,3:4); imagesc(cluster_image{level, ori});...
                    axis image; colormap(jet(256));
            end
            
        else
            %don't try and synthesise from a model that doesn't exist!
            warning(['Model ', syn_args.PathToTextureGMM, ' does not exist']); %#ok
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