function [synthesised_image cluster_image] = mb_gmm_tex_synthesis_goat(varargin)
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
          'UpperBand1',...
          'UpperBand2',...
          'UpperBand3',...
          'UpperBand4',...
          'UpperBand5',...
          'UpperBand6'...
          'MaxIdx',...
          },... % The optional arguments
          'LowerBand1', [], ...
          'LowerBand2', [], ...
          'LowerBand3', [], ...
          'LowerBand4', [], ...
          'LowerBand5', [], ...
          'LowerBand6', [], ...
          'PixelOrder', [],...
          'FilledImage', [], ...
          'ClusterImage', [], ...
          'WindowSize', 11,...
          'WindowSize2', 5,...
          'SynthesisMode', 'pixel-wise', ...
		  'MovieFilename', [], ...
		  'NormaliseClusterProbs', (1==1), ...
          'ForceCluster', 0, ...
          'SaveFrequency', 1000,...
          'SaveFile', []);

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

if isfield(model, 'WindowSize')
    args.WindowSize = model.WindowSize;
end
if isfield(model, 'WindowSize2')
    args.WindowSize2 = model.WindowSize2;
end

%Covert images to doubles for rest of algorithm
synthesised_image{1} = double(args.UpperBand1);
synthesised_image{2} = double(args.UpperBand2);
synthesised_image{3} = double(args.UpperBand3);
synthesised_image{4} = double(args.UpperBand4);
synthesised_image{5} = double(args.UpperBand5);
synthesised_image{6} = double(args.UpperBand6);
args.LowerBand1 = double(args.LowerBand1);
args.LowerBand2 = double(args.LowerBand2);
args.LowerBand3 = double(args.LowerBand3);
args.LowerBand4 = double(args.LowerBand4);
args.LowerBand5 = double(args.LowerBand5);
args.LowerBand6 = double(args.LowerBand6);

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
synthesised_image{1}(setdiff(args.PixelOrder, find(args.MaxIdx==1))) = NaN;
synthesised_image{2}(setdiff(args.PixelOrder, find(args.MaxIdx==2))) = NaN;
synthesised_image{3}(setdiff(args.PixelOrder, find(args.MaxIdx==3))) = NaN;
synthesised_image{4}(setdiff(args.PixelOrder, find(args.MaxIdx==4))) = NaN;
synthesised_image{5}(setdiff(args.PixelOrder, find(args.MaxIdx==5))) = NaN;
synthesised_image{6}(setdiff(args.PixelOrder, find(args.MaxIdx==6))) = NaN;

%if no cluster image supplied pre-allocate zeros - output image only shows
%assigned cluster in synthesised area
if isempty(args.ClusterImage)
    cluster_image = zeros(size(synthesised_image{1}));
    
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

% see whether we need to make the movie file
if ~isempty(args.MovieFilename)
    temp_im = [...
        synthesised_image{1} synthesised_image{2} synthesised_image{3};...
        synthesised_image{4} synthesised_image{5} synthesised_image{6}];
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
pixels_filled = 0;
tic; % Start a record of when we started.

[rows, cols] = ind2sub(size(synthesised_image{1}), args.PixelOrder);

% Here is the main algorithm.
%h = timebar(0,'Synthesising pixels for next onion skin...');
for ii = 1:pixels_to_fill
    
    row = rows(ii);
    col = cols(ii);
    
    upper_window1 = sample_window(synthesised_image{1}, args.WindowSize, row, col);
    upper_window2 = sample_window(synthesised_image{2}, args.WindowSize, row, col);
    upper_window3 = sample_window(synthesised_image{3}, args.WindowSize, row, col);
    upper_window4 = sample_window(synthesised_image{4}, args.WindowSize, row, col);
    upper_window5 = sample_window(synthesised_image{5}, args.WindowSize, row, col);
    upper_window6 = sample_window(synthesised_image{6}, args.WindowSize, row, col);
    
    if args.WindowSize2
        lower_window1 = sample_window(args.LowerBand1, args.WindowSize2, ceil(row/2), ceil(col/2));
        lower_window2 = sample_window(args.LowerBand2, args.WindowSize2, ceil(row/2), ceil(col/2));
        lower_window3 = sample_window(args.LowerBand3, args.WindowSize2, ceil(row/2), ceil(col/2));
        lower_window4 = sample_window(args.LowerBand4, args.WindowSize2, ceil(row/2), ceil(col/2));
        lower_window5 = sample_window(args.LowerBand5, args.WindowSize2, ceil(row/2), ceil(col/2));
        lower_window6 = sample_window(args.LowerBand6, args.WindowSize2, ceil(row/2), ceil(col/2));
    else
        lower_window1 = [];
        lower_window2 = [];
        lower_window3 = [];
        lower_window4 = [];
        lower_window5 = [];
        lower_window6 = [];
    end
    
    
    sampled_window = [upper_window1(:)' lower_window1(:)' ...
                      upper_window2(:)' lower_window2(:)' ...
                      upper_window3(:)' lower_window3(:)' ...
                      upper_window4(:)' lower_window4(:)' ...
                      upper_window5(:)' lower_window5(:)' ...
                      upper_window6(:)' lower_window6(:)'];    
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
    upper_window1 = sampled_window(1:args.WindowSize^2,1);
    upper_window2 = sampled_window(1:args.WindowSize^2,2);
    upper_window3 = sampled_window(1:args.WindowSize^2,3);
    upper_window4 = sampled_window(1:args.WindowSize^2,4);
    upper_window5 = sampled_window(1:args.WindowSize^2,5);
    upper_window6 = sampled_window(1:args.WindowSize^2,6);
    
    % Determine whether to place the whole window in the image, or
    % just the sampled pixel
    if pixel_wise_synthesis
        % Place the sampled pixel into the image.
        synthesised_image{1}(row, col) = upper_window1(centre_pixel_dim);
        synthesised_image{2}(row, col) = upper_window2(centre_pixel_dim);
        synthesised_image{3}(row, col) = upper_window3(centre_pixel_dim);
        synthesised_image{4}(row, col) = upper_window4(centre_pixel_dim);
        synthesised_image{5}(row, col) = upper_window5(centre_pixel_dim);
        synthesised_image{6}(row, col) = upper_window6(centre_pixel_dim);

    else
        % Patch-wise synthesis:
        % reshape top level of sampled_window into a square
        upper_window1 = reshape(upper_window1, [args.WindowSize args.WindowSize]);
        upper_window2 = reshape(upper_window2, [args.WindowSize args.WindowSize]);
        upper_window3 = reshape(upper_window3, [args.WindowSize args.WindowSize]);
        upper_window4 = reshape(upper_window4, [args.WindowSize args.WindowSize]);
        upper_window5 = reshape(upper_window5, [args.WindowSize args.WindowSize]);
        upper_window6 = reshape(upper_window6, [args.WindowSize args.WindowSize]);
        
        % now place the window in the image at the appropriate location
        [synthesised_image{1}] = ...
            mb_place_window(synthesised_image{1}, upper_window1, row, col);
        [synthesised_image{2}] = ...
            mb_place_window(synthesised_image{2}, upper_window2, row, col);
        [synthesised_image{3}] = ...
            mb_place_window(synthesised_image{3}, upper_window3, row, col);
        [synthesised_image{4}] = ...
            mb_place_window(synthesised_image{4}, upper_window4, row, col);
        [synthesised_image{5}] = ...
            mb_place_window(synthesised_image{5}, upper_window5, row, col);
        [synthesised_image{6}] = ...
            mb_place_window(synthesised_image{6}, upper_window6, row, col);
        
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
        temp_im = [...
            synthesised_image{1} synthesised_image{2} synthesised_image{3};...
            synthesised_image{4} synthesised_image{5} synthesised_image{6}];
        
        write_gif_from_colormap(temp_im, args.MovieFilename,...
            jet(256), 'append', movie_lims);
        clear temp_im;
        save_movie_percentages(1) = [];
    end
    
    % see if we need to save
    if ~(mod(pixels_filled, args.SaveFrequency)) && ~isempty(args.SaveFile);
        save(args.SaveFile, 'synthesised_image*');
    end
end

% save final image
if ~isempty(args.SaveFile)
    save(args.SaveFile, 'synthesised_image*');
end

% Report the generation of the synthetic image
disp('Progress: 100 percent complete.');
disp('A synthetic ergodic texture has been synthesised!');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%