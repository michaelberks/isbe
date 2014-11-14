function [synthesised_bands cluster_image] = mb_gmm_tex_synthesis_cls(varargin)
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
% 'TargetImage'
%   - 
%
% 'SampleImage2'
%   - the full path to an image (readable by Matlab's IMREAD) which
%   will be sampled to seed the synthetic texture (e.g. could be an
%   image from the training set used to build the GMM).
%
% 'FilledImage'
%   - a binary image that indicates where TargetImage has
%   been seeded (1's) and where the synthesis needs to operate
%   (zeroes). We use the sum of all the values in FilledImage to
%   determine if we need to automatically seed SYNTHESISED_IMAGE, so
%   if you want the seeding to be done automatically, you can just set
%   FilledImage to 0.
%
% Mandatory Arguments:
%
% 'PathToTextureGMM'
%   - the full path to a GMM model of texture, as built using the
%   BUILD_GM_TEXTURE_MODEL function.
%
% Optional Arguments:
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
          'UpperBands'...
          'ClsMap',...
          },... % The optional arguments
          'LowerBands', [], ...
          'PixelOrder', [],...
          'FilledImage', [], ...
          'ClusterImage', [], ...
          'WindowSize', 11,...
          'WindowSize2', 5,...
          'SynthesisMode', 'pixel', ...
		  'MovieFilename', [], ...
		  'NormaliseClusterProbs', (1==1), ...
          'ForceCluster', 0, ...
          'SaveFrequency', 1000,...
          'SaveFile', []);

% Use a proxy for the 'SynthesisMode' parameter, to avoid doing
% expensive string comparisons during run-time.
switch args.SynthesisMode
 case 'pixel'
  pixel_wise_synthesis_orig = (1==1);
 case 'patch'
  pixel_wise_synthesis_orig = (1==0);
 otherwise
  error(['The ' args.SynthesisMode ' mode is not supported.']);
end

% Load and check the models.
normal_model = u_load([args.PathToTextureGMM, '_cls_0']);
cls_model = u_load([args.PathToTextureGMM, '_cls_0']);

if args.NormaliseClusterProbs
  % Make sure that the ClusterProbs sum to unity
  normal_model.ClusterProbs = normal_model.ClusterProbs / sum(normal_model.ClusterProbs);
  cls_model.ClusterProbs = cls_model.ClusterProbs / sum(cls_model.ClusterProbs);
end
mb_gmm_check_model(normal_model);
mb_gmm_check_model(cls_model);

if isfield(normal_model, 'WindowSize')
    args.WindowSize = normal_model.WindowSize;
end
if isfield(normal_model, 'WindowSize2')
    args.WindowSize2 = normal_model.WindowSize2;
end

%see if we have an order inwhich to synthesis pixels, if not we have to
%make one from filled image
if isempty(args.PixelOrder)
    if ~isempty(args.FilledImage)
        cls_filled_image = args.FilledImage | ~args.ClsMap;
        normal_filled_image = args.FilledImage | args.ClsMap;
        
        cls_pixel_order = ...
            mb_compute_pixel_order(cls_filled_image, 'patch', args.WindowSize);
        normal_pixel_order = ...
            mb_compute_pixel_order(normal_filled_image, args.SynthesisMode, args.WindowSize);
        args.PixelOrder = [cls_pixel_order; normal_pixel_order]; 
    else
        error('Must supply either pixel order or filled image');
    end
end

%Discard the old values from region to synthesise in each sub-band;
upper_bands = cell(6,1);
for band = 1:6
    upper_bands{band} = args.UpperBands(:,:,band);
    upper_bands{band}(args.PixelOrder) = complex(NaN, NaN);
end

%Get size of each upper sub-band
[R C] = size(upper_bands{band});

%If we're using lower bands we need to interpolate these now
if args.WindowSize2
    
    %pre-allocate cell for new interpolated coarse bands
    lower_bands = cell(6,1);
    
    % Set up the expected phase shifts for interpolating each subband:
    w = [-3 -1; -3 -3; -1 -3; 1 -3; 3 -3; 3 -1]*pi/2.15;
    
    for band = 1:6
        temp = cpxinterp2(args.LowerBands(:,:,band), [-.25 .25], w(band,:),'spline');        
        lower_bands{band} = temp(1:R, 1:C);
    end
end

%if no cluster image supplied pre-allocate zeros - output image only shows
%assigned cluster in synthesised area
if isempty(args.ClusterImage)
    cluster_image = zeros(R, C);
    
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

% Find the dims of the window vectors that corresponds to the centre
% pixel.
centre_pixel_dim = .5*(args.WindowSize.^2+1);
dimensions = 2*(args.WindowSize^2 + args.WindowSize2^2);

% see whether we need to make the movie file
if ~isempty(args.MovieFilename)
    temp_im = [];
    for band = 1:6
        temp_im = [temp_im [real(upper_bands{band});imag(upper_bands{band})]]; %#ok
    end
    movie_lims(1) = min(temp_im(:));
    movie_lims(2) = max(temp_im(:));
    
    write_gif_from_colormap(temp_im, args.MovieFilename,...
            jet(256), 'overwrite', movie_lims);
    clear temp_im;
    save_movie_percentages = linspace(0, 1, 120);
end


%We do some manual progress feedback to the user, and set it up
%here. First, define a set of percentages to DISP at.
progress_report_percentages = 0.1:0.1:1;
pixels_to_fill = length(args.PixelOrder);
cls_pixels_to_fill = length(cls_pixel_order);
pixels_filled = 0;
tic; % Start a record of when we started.

[rows, cols] = ind2sub([R C], args.PixelOrder);

% Here is the main algorithm.
cls_count = 0;
normal_count = 0;
%h = timebar(0,'Synthesising pixels for next onion skin...');
for ii = 1:pixels_to_fill
    
    if ii <= cls_pixels_to_fill 
        %Use cls model
        model = cls_model;
        pixel_wise_synthesis = (1 == 0);
        cls_count = cls_count+1;
    else
        %Use normal model
        model = normal_model;
        pixel_wise_synthesis = pixel_wise_synthesis_orig;
        normal_count = normal_count+1;
    end
    
    row = rows(ii);
    col = cols(ii);
    sampled_window = zeros(1, 6*dimensions);
    for band = 1:6
        sample_win_re1 = ...
            sample_window(real(upper_bands{band}), args.WindowSize, row, col);
        sample_win_im1 = ...
            sample_window(imag(upper_bands{band}), args.WindowSize, row, col);
    
        if args.WindowSize2
            sample_win_re2 = ...
                sample_window(real(lower_bands{band}), args.WindowSize2, row, col);
            sample_win_im2 = ...
                sample_window(imag(lower_bands{band}), args.WindowSize2, row, col);
        else
            sample_win_re2 = [];
            sample_win_im2 = [];
        end
    
    
        sampled_window((band-1)*dimensions + 1:band*dimensions) =...
                [sample_win_re1(:)' sample_win_im1(:)' sample_win_re2(:)' sample_win_im2(:)'];
    end
    
    % Sample from the model:
    
    %First condition the model using known dimensions
    conditioned_model =...
        mb_gmm_condition(model, sampled_window);
    
    %See if we've been told to choose a specific cluster from the model
    force_cluster = cluster_image(row, col);
    if isnan(force_cluster);
         force_cluster = [];
        %Just test to see if choosing the most likely cluster works better
        %[dummy force_cluster] = max(conditioned_model.ClusterProbs);
    end
    
    % now sample from the conditioned model
    [p chosen_cluster] = mb_gmm_sample('Model', conditioned_model, 'ForceCluster', force_cluster);
    cluster_image(row, col) = chosen_cluster;
   
    %sometimes returns zero complex part - not sure why but scrap it
    p = real(p);
    
    % now we have p, place it in the sampled window in place of the NaN's
    % unstandardising and re-complexing as we do
    sampled_window(isnan(sampled_window)) = p;
    sampled_window = reshape(sampled_window, [], 6);
    
    for band = 1:6
        sample_win_re1 = sampled_window(1:args.WindowSize^2, band);
        sample_win_im1 = sampled_window(args.WindowSize^2 + 1:2*args.WindowSize^2, band);
        
    
        % Determine whether to place the whole window in the image, or
        % just the sampled pixel
        if pixel_wise_synthesis
            % Place the sampled pixel into the image.
            upper_bands{band}(row, col) =...
                complex(sample_win_re1(centre_pixel_dim), sample_win_im1(centre_pixel_dim));

        else
            % Patch-wise synthesis:
            % reshape top level of sampled_window into a square
            sample_win_re1 = reshape(sample_win_re1, [args.WindowSize args.WindowSize]);
            sample_win_im1 = reshape(sample_win_im1, [args.WindowSize args.WindowSize]);

            % now place the window in the image at the appropriate location
            [upper_bands{band}] = ...
                mb_place_window(upper_bands{band}, complex(sample_win_re1, sample_win_im1), row, col);
        
        end
    end
    
    % Update the record of pixel filled
    pixels_filled = pixels_filled +  1;
    
    % Report progress
    progress_so_far = pixels_filled / pixels_to_fill;
    if progress_so_far > progress_report_percentages(1)
      disp(['Progress: ' num2str(floor(100 * progress_so_far)) ' percent' ...
		    ' complete.']);
      disp(['--Time: ' datestr(now)]);
      % Remove the first report point from the progress_report_percentages
      progress_report_percentages(1) = [];
    end
    
    % see if we have to update the movie
    if ~isempty(args.MovieFilename) && progress_so_far >= save_movie_percentages(1);
        temp_im = [];
        for band = 1:6
            temp_im = [temp_im [real(upper_bands{band});imag(upper_bands{band})]]; %#ok
        end
        write_gif_from_colormap(temp_im, args.MovieFilename,...
            jet(256), 'append', movie_lims);
        clear temp_im
        save_movie_percentages(1) = [];
    end
    
    % see if we need to save
    if ~(mod(pixels_filled, args.SaveFrequency)) && ~isempty(args.SaveFile);
        save(args.SaveFile, 'upper_bands*');
    end
end

% save final movie frames
if ~isempty(args.MovieFilename)
    temp_im = [];
    for band = 1:6
        temp_im = [temp_im [real(upper_bands{band});imag(upper_bands{band})]]; %#ok
    end
    write_gif_from_colormap(temp_im, args.MovieFilename,...
        jet(256), 'append', movie_lims);
    clear temp_im
end

%Reshape synthesised bands into array structure
synthesised_bands = zeros(R, C, 6);
for band = 1:6
    synthesised_bands(:,:,band) = upper_bands{band};
end

% save final image
if ~isempty(args.SaveFile)
    save(args.SaveFile, 'synthesised_bands*');
end

% Report the generation of the synthetic image
disp('Progress: 100 percent complete.');
disp('A synthetic ergodic texture has been synthesised!');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%