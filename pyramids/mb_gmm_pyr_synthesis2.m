function [synthesised_image, pyramid] = mb_gmm_pyr_synthesis2(varargin)
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
          'PathToPyramidGMM', ...
          'SampleImage',...
          'FilledImage'...
          },... % The optional arguments
          'SynthesisMode', 'simple', ...
		  'MovieFilename', [], ...
		  'NormaliseClusterProbs', (1==1), ...
		  'SaveFrequency', 1000,...
          'SaveFile', [mberksroot, 'background\syn\temp\syn_image'],...
          'PyrLevels', 5,...
          'PyrOrientations', 5,...
          'CutOffLevel', 2);

% Use a proxy for the 'SynthesisMode' parameter, to avoid doing
% expensive string comparisons during run-time.
switch args.SynthesisMode
 case 'simple'
  
 case 'advanced'
  
 otherwise
  error(['The ' args.SynthesisMode ' mode is not supported.']);
end

% Load and check the model.
model = u_load(args.PathToPyramidGMM);

if args.NormaliseClusterProbs
  % Make sure that the ClusterProbs sum to unity
  model.ClusterProbs = model.ClusterProbs / sum(model.ClusterProbs);
end
mb_gmm_check_model(model);


% We need SampleImage and args.SeededImage to be uint8's because we're assuming pixel values in the
% range 0-255.
if ~isa(args.SampleImage, 'uint8')
    args.SampleImage = uint8(args.SampleImage);
    warning('Automatically converted the sample image (SampleImage) to a uint8') %#ok
end
% convert SampleImage and args.SeededImage to doubles, maintain the pixel value ranges (0-255)
args.SampleImage = double(args.SampleImage); % SampleImage is now a double

% see whether we need to make the movie file
if ~isempty(args.MovieFilename)
    movie = avifile(args.MovieFilename, 'compression', 'none'); %#ok
end

% get row/col subscripts of unfilled pixels
[p_rows p_cols] = find(~args.FilledImage);
%p_idx = sub2ind(size(args.SampleImage), p_rows, p_cols);

% Calculate the pyramid for the image to be synthesised
[pyramid p_sizes] = buildSFpyr(args.SampleImage, args.PyrLevels, args.PyrOrientations-1); %

% Convert pyrmaid into Rose form
pyramid = change_pyr_structure(pyramid, args.PyrLevels, args.PyrOrientations, p_sizes);
    
% Get pyramid coefficients for unfilled pixels
pyr_params = mb_get_pyramid_coefficients(pyramid, [p_rows, p_cols]);

% Discard finest levels of pyramid by replacing params with Nan
cut_off = args.CutOffLevel*args.PyrOrientations + 1;
pyr_params(:, 2:cut_off) = NaN; 

tic; % Start a record of when we start actually working on the pixels
%We do some manual progress feedback to the user, and set it up
%here. First, define a set of percentages to DISP at.
% progress_report_percentages = 0.1:0.1:1;
% pixels_to_fill = sum(~args.FilledImage(:));

for level = args.CutOffLevel+1:-1:2
    
    % Calculate indices for that level and make unique
    new_rows = ceil(p_rows/2^(level-2));
    new_cols = ceil(p_cols/2^(level-2));
    new_idx = sub2ind(size(pyramid{level,1}), new_rows, new_cols);
    
    [new_idx orig_idx something_else] = unique(new_idx); 
    %new_subs is 2-column matrix of unique indices, new_idx is column
    %vector specifying the original indices of the new subs
    
    % calculate the start and end dimensions of the orientation parameters
    % at this level
    start_dim = 2 + (level-2)*args.PyrOrientations;
    end_dim = 1 + (level-1)*args.PyrOrientations;
    
    level_pixels = length(new_idx);
    level_params = zeros(level_pixels, sum(isnan(pyr_params(1,:))));
    
    % loop through each pixel conditionally sampling new pyramid parameters  
    for jj = 1:level_pixels

        % Condition model based on known components  
        conditioned_model = mb_gmm_condition(model, pyr_params(orig_idx(jj),:));

        % now sample from the conditioned model
        level_params(jj,:) = mb_gmm_sample('Model', conditioned_model); %Horrible, horrible hack!!!
        
    end
    %place p back in pyr_params. p is full params vector, we only want 
    %the orientations params at this level 
    pyr_params(:,start_dim:end_dim) = level_params(something_else,end-args.PyrOrientations+1:end);

    %Also put the params in the correct level of the full pyramid
    for ori = 1:args.PyrOrientations
        pyramid{level, ori}(new_idx) = level_params(:,ori);
    end

end

% put pyrimad back into Simonecelli form
pyramid = change_pyr_structure(pyramid, args.PyrLevels, args.PyrOrientations);

% reconstruct image from pyrimad
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
save(args.SaveFile, 'synthesised_image');

% Report the generation of the synthetic image
disp('Progress: 100 percent complete.');
disp('A synthetic ergodic texture has been synthesised!');

%
end %main function end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function new_pyr = change_pyr_structure(pyramid, num_levels, num_oris, pind)
    
if iscell(pyramid)
    %Pyrmaid is in Rose form - change to Simmoncelli
    %new pyramind is single column vector; could pre-allocate but can't be
    %arsed - we only do this once so performance isn't key;
    new_pyr = [];
    for lev = 1:num_levels+2
        if lev == 1 || lev == num_levels+2
            % high/low band pass - no sub-band orientations
            new_pyr = [new_pyr; pyramid{lev,1}(:)]; %#ok
        
        else
            %levels comprised of num_oris sub-band orientations
            for ori = 1:num_oris
                new_pyr = [new_pyr; pyramid{lev, ori}(:)]; %#ok
            end
        end
    end
        
    
else
    %Pyramid is in Simmoncelli form - change to Rose
    new_pyr = cell(num_levels+2, num_oris);
    band_idx = cumsum(prod(pind, 2));
    
    for lev = 1:num_levels+2
        if lev == 1
            % high band pass - no sub-band orientations
            new_pyr{lev, 1} = reshape(pyramid(1:band_idx(1)), pind(1,1), pind(1,2));
        elseif lev == num_levels+2
            %low band pass - no sub-band orientations
            new_pyr{lev, 1} = reshape(pyramid(band_idx(end-1)+1:end), pind(end,1), pind(end,2));
        else
            %levels comprised of num_oris sub-band orientations
            for ori = 1:num_oris
                band = num_oris*(lev-2) + ori + 1;
                new_pyr{lev, ori} = reshape(pyramid(band_idx(band-1)+1:band_idx(band)),...
                    pind(band,1), pind(band,2));
            end
        end
    end
end

end % end of sub-function
















