function [synthesised_image, pyramid, cluster_image] = mb_gmm_pyr_synthesis(varargin)
%
%MB_GMM_PYR_SYNTHESIS Synthesise an image using an ergodic GMM.
%
% 
%
% MB_GMM_PYR_SYNTHESIS uses the U_PACKARGS interface function
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
%   MB_GMM_PYR_SYNTHESIS returns SYNTHESISED_IMAGE
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
          'TargetPyramid', [],...
          'ClusterImage', [], ...
          'SynthesisMode', 'simple', ...
		  'SaveFile', [],...
          'PyrLevels', 5,...
          'PyrOrientations', 5,...
          'CutOffLevel', 5,...
          'ForceCluster', 0, ...
          'Plot', 0,...
          'ConditionLevels', 1,...
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
    
    % Calculate the pyramid for the image to be synthesised
    [pyramid p_sizes] = mb_buildSFpyr(args.TargetImage, args.PyrLevels, args.PyrOrientations-1); %
    
    % Convert pyramid into Rose form
    pyramid = mb_change_pyramid_form(pyramid, p_sizes);

elseif ~isempty(args.TargetPyramid)
    pyramid = args.TargetPyramid;
    args = rmfield(args, 'TargetPyramid');
else
    error('Either a target image or target pyramid must be supplied');
end

% get row/col subscripts of unfilled pixels
[p_rows p_cols] = find(~args.FilledImage);
%p_idx = sub2ind(size(args.TargetImage), p_rows, p_cols);

tic; % Start a record of when we start actually working on the synthesis

%pre-allocate cluster_image for output
cluster_image = cell(args.PyrLevels+2, args.PyrOrientations);

for level = args.CutOffLevel:-1:2
    
    if level > 1
        % Calculate indices for that level and make unique
        new_rows = ceil(p_rows/2^(level-2));
        new_cols = ceil(p_cols/2^(level-2));
        new_idx = sub2ind(size(pyramid{level,1}), new_rows, new_cols);

        %Don't need to do unique now
        %[new_idx orig_idx something_else] = unique(new_idx);
        %new_subs is 2-column matrix of unique indices, new_idx is column
        %vector specifying the original indices of the new subs

        %make filled image from new_idx
        filled_image = ones(size(pyramid{level,1}));
        filled_image(new_idx) = 0;
    else
        filled_image = args.FilledImage;
    end
    
    syn_args.PixelOrder = mb_compute_pixel_order(filled_image);
    
    %make sure uint8 mode is turned off since we're synthesising doubles
    syn_args.Uint8Mode = 0;
    
    %set force clusters flag
    syn_args.ForceCluster = args.ForceCluster;
    
    for ori = 1:args.PyrOrientations
        
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
            syn_args.TargetImage = pyramid{level, ori};
            
            if ~isempty(args.ClusterImage)
                syn_args.ClusterImage = args.ClusterImage{level, ori};
            end
            
            %Now synthesise the new texture for the band and save back to
            %pyramid
            %Work out whether we're conditioning on the lower level
            if args.ConditionLevels
                %Use mb_gmm_tex_synthesis2
                %Need to give lower pyramid level as sample image
                if level + 1 == size(pyramid, 1)
                    syn_args.SampleImage2 = pyramid{level+1, 1};
                else
                    syn_args.SampleImage2 = pyramid{level+1, ori};
                end
                
                [pyramid{level, ori} cluster_image{level, ori}] = ...
                        mb_gmm_tex_synthesis2(syn_args);
            else
                %Use mb_gmm_tex_synthesis
                [pyramid{level, ori} cluster_image{level, ori}] = ...
                    mb_gmm_tex_synthesis(syn_args);
            end

            %display new texture if required
            if args.Plot
                figure; imagesc(pyramid{level, ori});...
                    axis image; colormap(jet(256));
                figure; imagesc(cluster_image{level, ori});...
                    axis image; colormap(jet(256));
            end
            
        else
            %don't try and synthesise from a model that doesn't exist!
            warning(['Model ', syn_args.PathToTextureGMM, ' does not exist']); %#ok
        end
        if level == 1
            break;
        end
        
    end

end

% put pyrimad back into Simonecelli form
[pyramid p_sizes] = mb_change_pyramid_form(pyramid);

% reconstruct image from pyramid
synthesised_image = uint8(reconSFpyr(pyramid, p_sizes));

%but actually we prefer Rose's form to save it in!
pyramid = mb_change_pyramid_form(pyramid, p_sizes);

filled_image = args.FilledImage; %#ok

% save final image
if ~isempty(args.SaveFile)
    save(args.SaveFile, 'synthesised_image', 'pyramid', 'cluster_image', 'filled_image');
end

% Report the generation of the synthetic image
disp('Progress: 100 percent complete.');
disp('A modified mamogram region has been synthesised!');

%Main function end