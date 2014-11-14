function [synthesised_image, pyramid] = mb_gmm_pyr_synthesis3(varargin)
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
% 'SeededImage'
%   - an image that is seeded and will be populated with
%   likely pixels to form synthesised_image. If seeding is to be
%   performed automatically, then set SeededImage to zeros(sz), where sz
%   is a vector that specifies the desired size of SYNTHESISED_IMAGE,
%   and set FilledImage to zero. Defaults to zeros(200,200) --
%   i.e. a 200x200 automatically-seeded image
%
% 'SampleImage'
%   - the full path to an image (readable by Matlab's IMREAD) which
%   will be sampled to seed the synthetic texture (e.g. could be an
%   image from the training set used to build the GMM).
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
% 'SeedPlacement'
%   - specified where the seed is placed. Can be one of:
%    + 'centre' (Default)
%    + 'bottom-right'
%    + 'bottem-left'
%    + 'top-right'
%    + 'top-left'
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
args = u_packargs(varargin, '0', ... % strict mode
		  {...  % The mandatory arguments
          'PathToPyramidGM', ...
          'FilledImage'...
          },... % The optional arguments
          'SampleImage', [],...
          'SamplePyramid', [],...
          'SynthesisMode', 'patch-wise', ...
		  'MovieFilename', [], ...
		  'SaveFrequency', 1000,...
          'SaveFile', [],...
          'PyrLevels', 5,...
          'PyrOrientations', 5,...
          'CutOffLevel', 4,...
          'WindowSize', 15);

      clear varargin

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

%Load the model
model = u_load(args.PathToPyramidGM);

if ~isempty(args.SampleImage)
    % convert SampleImage to doubles, maintain the pixel value ranges (0-255)
    args.SampleImage = double(args.SampleImage); % SampleImage is now a double
    
    % Calculate the pyramid for the image to be synthesised
    [pyramid p_sizes] = buildSFpyr(args.SampleImage, args.PyrLevels, args.PyrOrientations-1); %
    
    % Convert pyramid into Rose form
    pyramid = mb_change_pyramid_form(pyramid, p_sizes);

elseif ~isempty(args.SamplePyramid)
    pyramid = args.SamplePyramid;
    args = rmfield(args, 'SamplePyramid');
else
    error('Either a sample image or sample pyramid must be supplied');
end


% see whether we need to make the movie file
if ~isempty(args.MovieFilename)
    movie = avifile(args.MovieFilename, 'compression', 'none'); %#ok
end

% get row/col subscripts of unfilled pixels
[p_rows p_cols] = find(~args.FilledImage);
%p_idx = sub2ind(size(args.SampleImage), p_rows, p_cols); 

centre_pixel_dim = .5*(args.WindowSize.^2+1);

tic; % Start a record of when we start actually working on the pixels
%We do some manual progress feedback to the user, and set it up
%here. First, define a set of percentages to DISP at.
% progress_report_percentages = 0.1:0.1:1;
% pixels_to_fill = sum(~args.FilledImage(:));

for level = args.CutOffLevel:-1:2
    
    % Calculate indices for that level
    new_rows = ceil(p_rows/2^(level-2));
    new_cols = ceil(p_cols/2^(level-2));
    new_idx = sub2ind(size(pyramid{level,1}), new_rows, new_cols);
    
    %use to make new mask
    mask = ones(size(pyramid{level, 1}));
    mask(new_idx) = 0;
    
    %Discard parameters from pyramid
    for ori = 1:args.PyrOrientations
        pyramid{level, ori}(new_idx) = NaN;
    end
    
    % get a list of the unfilled pixels
    unfilled_list = get_unfilled_perimeter_pixels(mask);
    
    %Compute cutoff point for coarse levels of pyramid
    cut_off = (level-1)*args.PyrOrientations + 2;
    
    % loop through each pixel conditionally sampling new pyramid parameters  
    while ~isempty(unfilled_list)
        
        %get rows/cols from indices
        [rows cols] = ind2sub(size(pyramid{level,1}), unfilled_list);
        
        %Get pyramid coeffs for rows/cols (take note of level!)
        pyr_params = mb_get_pyramid_coefficients(pyramid, 2^(level-2)*[rows, cols]);
        
        %compute random permutation indices
        rand_idx = randperm(length(unfilled_list));
        
        for jj = 1:length(unfilled_list)
            
            % get row and col from random list
            row = rows(rand_idx(jj));
            col = cols(rand_idx(jj));
            
            %pre-allocate sample window
            sampled_window = zeros(1, args.PyrOrientations*args.WindowSize^2);
            
            %take a sample window from pixel in each orientated sub-band
            for ori = 1:args.PyrOrientations
                sampled_window((ori-1)*args.WindowSize^2 + 1:ori*args.WindowSize^2) = ...
                    reshape(sample_window(pyramid{level, ori}, ...
                          args.WindowSize, row, col), 1, []);    
            end
            
            %add coarse level pyr_params to the data
            sampled_window = [pyr_params(rand_idx(jj), cut_off:end),...
                sampled_window]; %#ok
            
            %Standardise data
            [sampled_window] = st_standardise_data(sampled_window,...
                model.St_Means{level}, model.St_Stds{level});
            
            % Condition model based on known components
            [cond_mean cond_covar] = ...
                condition_gaussian(model.Means{level}, model.Covars{level}, sampled_window);

            %sample new coefficients using conditioned model
            new_coeffs = real(sample_from_normal(cond_mean, cond_covar, 1));

            %Now try and put everything back in the right place!!
            sampled_window(isnan(sampled_window)) = new_coeffs;
            
            %Unstandardise data
            [sampled_window] = st_unstandardise_data(sampled_window,...
                model.St_Means{level}, model.St_Stds{level});
            
            %get rid pyr_params
            sampled_window = ...
                sampled_window(end - args.PyrOrientations*args.WindowSize^2 + 1: end);
            
            for ori = 1:args.PyrOrientations
                % reshape sampled_window into a square
                single_window = reshape(...
                    sampled_window((ori-1)*args.WindowSize^2 + 1:ori*args.WindowSize^2),...
                    [args.WindowSize args.WindowSize]);

                % Determine whether to place the whole window in the image, or
                % just the sampled pixel (i.e. decide between patch- and
                % pixel-wise synthesis).
                if pixel_wise_synthesis
                  % pixel-wise synthesis:
                  % Place the sampled pixel into the image.
                  pyramid{level, ori}(row, col) = single_window(centre_pixel_dim);
                  % Update the record of where we've been.
                  mask(row, col) = 1;
                else
                  % Patch-wise synthesis:
                  % now place the window in the image at the appropriate location
                  [pyramid{level, ori}] =...
                      mb_place_window(pyramid{level, ori}, single_window, row, col);
                  %pixels_filled = pixels_filled +  length(p);                  
                end
                pixels_filled = pixels_filled +  1;
            end
        end
        % get new list of the unfilled pixels
        unfilled_list = get_unfilled_perimeter_pixels(mask);
    end
end

% put pyrimad back into Simonecelli form
[pyramid p_sizes] = mb_change_pyramid_form(pyramid);

% reconstruct image from pyramid
synthesised_image = uint8(reconSFpyr(pyramid, p_sizes));

% % Report progress
% progress_so_far = jj / pixels_to_fill;
% if progress_so_far > progress_report_percentages(1)
%   disp(['Progress: ' num2str(floor(100 * progress_so_far)) ' percent' ...
%         ' complete.']);
%   disp(['--Time: ' datestr(now)]);
%   % Remove the first report point from the progress_report_percentages
%   progress_report_percentages(1) = [];
% end

% save final image
if ~isempty(args.SaveFile)
    save(args.SaveFile, 'synthesised_image', 'pyramid');
end

% Report the generation of the synthetic image
disp('Progress: 100 percent complete.');
disp('A synthetic ergodic texture has been synthesised!');

%
end %main function end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [unfilled_list] = get_unfilled_perimeter_pixels(filled_image)
%
% Get a list of unfilled pixels neighbouring the pixels already filled.
% This is returned as a list of indices -- use ind2sub to get the image
% locations back.

unfilled_image = imdilate(filled_image, strel('disk', 1));
unfilled_image = double(unfilled_image) - double(filled_image);
unfilled_list = find(unfilled_image==1); 
end















