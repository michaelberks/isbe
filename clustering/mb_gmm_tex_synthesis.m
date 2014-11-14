function [synthesised_image cluster_image] = mb_gmm_tex_synthesis(varargin)
%
%MB_GMM_TEX_SYNTHESIS Synthesise an image using an ergodic GMM.
%
% This function synthesises an ergodic texture using a Gaussian
% Mixture Model of texture. The code can operate in pixel-wise
% or patch-wise modes. The order in which pixels/patches are
% populated uses the 'onion-skin' ordering -- the 'fill-in' scheme
% does not work very well.
%
% MB_GMM_TEX_SYNTHESIS uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%
% 'PathToTextureGMM'
%   - the full path to a GMM model of texture, as built using the
%   BUILD_GM_TEXTURE_MODEL function.
%
% Optional Arguments:
%
% 'TargetImage'
%   - 
%
% 'FilledImage'
%   - a binary image that indicates where SeededImage has
%   been seeded (1's) and where the synthesis needs to operate
%   (zeroes). We use the sum of all the values in FilledImage to
%   determine if we need to automatically seed SYNTHESISED_IMAGE, so
%   if you want the seeding to be done automatically, you can just set
%   FilledImage to 0.
%
% 'SynthesisMode'
%   - specifies whether synthesis is performed on a pixel- or
%   patch-wise basis. Use 'pixel-wise' or 'patch-wise'. Defaults to
%   'pixel-wise'.
%
% 'MovieFilename' 
%   - this function can capture a movie of the synthesis
%   process. Specify a filename to save the movie under to enable this.
%
% 'NormaliseClusterProbs'
%   - specifies whether to force GMM component probabilities to sum to
%   unity. Set to true to force normalisation (the default), or
%   false. (In early versions of this code, components of the GMM were
%   not assigned probabilities, but counts of the number of training
%   points assigned. To ensure that we have a reasonably valid GMM, we
%   perform a check of it using MB_GMM_CHECK_MODEL. Apart from other
%   things, this checks that the model's component probabilities sum
%   to unity. If this check fails, as it would do for early models,
%   then an error is thrown. This parameter can be used to ensure that
%   this does not happen, or to be paranoid and make sure that the GMM
%   is as valid as we can.)
%
%
% Return Value:
%
%   MB_GMM_TEX_SYNTHESIS returns SYNTHESISED_IMAGE, a synthetic texture.
%
% Example Usage:
%
%   % Generate a synthetic texture of default size.
%   syn_texture = mb_gmm_tex_synthesis(...
%                 'PathToTextureGMM', '~/my-model.mat', ...
%                 'PathToSampleImage', '~/my-seed-image.bmp');
%
% Known issues:
%
%   GMMs built in the natural space and where synthesis is performed
%   pixel-wise produce the best synthetic textures. Patch-wise
%   textures can be synthetised faster, but there is greater
%   opportunity for synthesis to go wrong. GMMs built in PCA spaces
%   can also produce poor textures; apart from having more confidence
%   in the density estimation task when building the model, there is
%   not real justification for using PCA GMMs.
%
% References:
%
%   This code is the implementation of the work described in
%   C. J. Rose and C. J. Taylor. A Statistical Model of Texture for
%   Medical Image Synthesis and Analysis. In Proc. Medical Image
%   Understanding and Analysis, pages 1-4, 2003.
%
% See also:
%
%   BUILD_GM_TEXTURE_MODEL, CJR_ANALYSE_TEXTURE


% Unpack the arguments:
args = u_packargs(varargin, 'strict', ... % strict mode
		  {...  % The mandatory arguments
          'PathToTextureGMM', ...
          'TargetImage'...
          },... % The optional arguments
          'PixelOrder', [],...
          'FilledImage', 0, ...
          'ClusterImage', [], ...
          'SynthesisMode', 'pixel-wise', ...
		  'MovieFilename', [], ...
		  'NormaliseClusterProbs', (1==1), ...
          'ForceCluster', 0, ...
		  'SaveFrequency', 1000,...
          'SaveFile', [],...
          'Uint8Mode', 0,...
          'ComplexMode', 0);

% Use a proxy for the 'SynthesisMode' parameter, to avoid doing
% expensive string comparisons during run-time.
switch args.SynthesisMode
 case 'pixel-wise'
  pixel_wise_synthesis = (1==1);
 case 'patch-wise'
  pixel_wise_synthesis = (1==0);
 otherwise
  error(['The ' args.SynthesisMode ' mode is not supported.']);
end

% Load and check the model.
model = u_load(args.PathToTextureGMM);
if args.NormaliseClusterProbs
  % Make sure that the ClusterProbs sum to unity
  model.ClusterProbs = model.ClusterProbs / sum(model.ClusterProbs);
end
mb_gmm_check_model(model);

%If we're using uint8 mode make sure images are in rang 0-255
if args.Uint8Mode
    args.TargetImage = uint8(args.TargetImage);
    warning('Automatically converted the seed image (args.TargetImage) to a uint8') %#ok

end

%Covert images to doubles for rest of algorithm
synthesised_image = double(args.TargetImage);

%see if we have an order inwhich to synthesis pixels, if not we have to
%make one from filled image
if isempty(args.PixelOrder)
    if ~isempty(args.FilledImage)
        args.PixelOrder = mb_compute_pixel_order(args.FilledImage);
    else
        error('Must supply either pixel order or filled image');
    end
end

%Discard the old values from region to synthesise;
synthesised_image(args.PixelOrder) = NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PCA stuff - ignore for now
% Determine if we are using a model built in a PCA space.
% pca_space_model = (1==0); % set to false, assume we are not dealing with a PCA model
% if isfield(model, 'pca_info_struct')
% 	% the field exists, but is the model in the PCA space or the natural data space?
% 	if length(model.Means(1,:)) == size(model.pca_info_struct.Eigenvectors,2)
% 		% the model is in the PCA space
% 		pca_space_model = (1==1); % set to true
% 	end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if no cluster image supplied pre-allocate zeros - output image only shows
%assigned cluster in synthesised area
if isempty(args.ClusterImage)
    cluster_image = zeros(size(synthesised_image));
    
    % Can't force clusters with no cluster image supplied
    if args.ForceCluster
        warning('Can''t force clusters with no cluster image.'); %#ok
        args.ForceCluster = 0;
    end
else
    %use cluster image given
    cluster_image = args.ClusterImage;
end

%if we want to choose cluster then set cluster values to NaN in synthesised
%region
if ~args.ForceCluster
    cluster_image(args.PixelOrder) = NaN;
end

% Find the dim of the window vectors that corresponds to the centre
% pixel.
window_size = sqrt(length(model.Means(1,:)));
centre_pixel_dim = .5*(window_size.^2+1);

% see whether we need to make the movie file
if ~isempty(args.MovieFilename)
    write_gif_from_colormap(synthesised_image, args.MovieFilename,...
            jet(256), [-4 4], 'overwrite');
    save_movie_percentages = linspace(0, 1, 120);
end


%We do some manual progress feedback to the user, and set it up
%here. First, define a set of percentages to DISP at.
progress_report_percentages = 0.1:0.1:1;
pixels_to_fill = length(args.PixelOrder);
pixels_filled = 0;
tic; % Start a record of when we started.

[rows, cols] = ind2sub(size(synthesised_image), args.PixelOrder);

% Here is the main algorithm.
%h = timebar(0,'Synthesising pixels for next onion skin...');
for ii = 1:pixels_to_fill
    
    row = rows(ii);
    col = cols(ii);
    
    sampled_window = reshape(sample_window(synthesised_image, ...
					  window_size, row, col), 1, []);
                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PCA stuff - ignore for now
%     % condition on what we do know
%     if pca_space_model
%         % the model is in the PCA space
%         model.PCA_GMMConditioningFunc = @pca_gm_tex_model_gmm_cond; 
%         % create the function handle that will perform conditioning on PCA GMMs
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Sample from the model:
    % Condition model using known pixels
    conditioned_model = mb_gmm_condition(model, sampled_window);
    
    %See if we've been told to choose a specific cluster from the model
    force_cluster = cluster_image(row, col);
    if isnan(force_cluster);
        force_cluster = [];
%         %Just test to see if choosing the most likely cluster works better
%         [dummy force_cluster] = max(conditioned_model.ClusterProbs);
    end
    
    % now sample from the conditioned model
    [p chosen_cluster] = mb_gmm_sample('Model', conditioned_model, 'ForceCluster', force_cluster);
    cluster_image(row, col) = chosen_cluster;
    
    if args.Uint8Mode
        % correct for the possibility that p may be a real number, and we want integers in the range 0-255
        p = double(uint8(p));
    elseif ~args.ComplexMode
        %sometimes returns zero complex part - not sure why but scrap it
        p = real(p);
    end
    
    % now we have p, place it in the sampled window in place of the NaN's
    sampled_window(isnan(sampled_window)) = p;
    
    % reshape sampled_window into a square
    sampled_window = reshape(sampled_window, [window_size window_size]);
    
    % Determine whether to place the whole window in the image, or
    % just the sampled pixel
    if pixel_wise_synthesis
      % Place the sampled pixel into the image.
      synthesised_image(row, col) = sampled_window(centre_pixel_dim);        
    else
      % Patch-wise synthesis:
      % now place the window in the image at the appropriate location
      [synthesised_image] = ...
          mb_place_window(synthesised_image, sampled_window, row, col); 
    end

    % Update the record of where we've been.
    pixels_filled = pixels_filled +  1;

%     % see if we have to update the movie
%     if ~isempty(args.MovieFilename)
%         frame(:,:,1) = uint8(synthesised_image);
%         frame(:,:,2) = uint8(synthesised_image);
%         frame(:,:,3) = uint8(synthesised_image); % this makes an NxMx3 uint8 frame
%         movie = addframe(movie, frame);
       
%     end
    
    % Report progress
    progress_so_far = pixels_filled / pixels_to_fill;
    
    if progress_so_far >= progress_report_percentages(1)
      disp(['Progress: ' num2str(floor(100 * progress_so_far)) ' percent' ...
		    ' complete.']);
      disp(['--Time: ' datestr(now)]);
      % Remove the first report point from the progress_report_percentages
      progress_report_percentages(1) = [];
    end
    
    % see if we have to update the movie
    if ~isempty(args.MovieFilename) && progress_so_far >= save_movie_percentages(1);
        write_gif_from_colormap(synthesised_image, args.MovieFilename,...
            jet(256), [-4 4], 'append');
        save_movie_percentages(1) = [];
    end
    
    % see if we need to save
    if ~(mod(pixels_filled, args.SaveFrequency)) && ~isempty(args.SaveFile)
        save(args.SaveFile, 'synthesised_image');
    end
end
% close(h);

% see if we have to close the movie fie
% if ~isempty(args.MovieFilename)
%     movie = close(movie); %#ok
% end

% save final image
if ~isempty(args.SaveFile)
    save(args.SaveFile, 'synthesised_image');
end

% Report the generation of the synthetic image
%disp('Progress: 100 percent complete.');
disp('A synthetic ergodic texture has been synthesised!');