function [synthesised_image, pixel_list, similarity_image] = efros_tex_syn(varargin)
%
% EFROS_TEX_SYN   CJR's implementation of the Efros and Leung algorithm.
%                   adapted by MAB
% EFROS_TEX_SYN uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%
% 'SampleImage'
%   - A sample of the texture that you want to synthesise.
%
% 'SeedImage'
%   - An image that will form the basis of the completed
%   texture. This image may already have some texture content, or
%   it can be completely blank. If you want this image to be seeded
%   automatically, then set all pixels to zero and set the
%   FilledImage parameter to zero.
%
% 'FilledImage' 
%   - A binary image that indicated where SeedImage
%   has been seeded (1's) and where the synthesis needs to
%   operate. Also read the info on SeedImage, above.
%
% 'WindowSize'
%   - the length of one side of the window. This must be an odd
%   number. For example, try using a value of 11.
%
%
% Optional Arguments:
%
% 'Debug'
%   - set this to 1 to perform args.Debugging.
%
% Return Value:
%
% EFROS_TEX_SYN returns a synthetic texture, with pixels
% populated by the Efros and Leung algorithm.
%
% Example Usage:
%
%   syn_image = cjr_efros_tex_syn('SampleImage', sample, ...
%      'SeedImage', zeros(30), FilledImage, 0);
%
% Known issues:
%
%  None so far.
%
% References:
%
%  A. A. Efros and T. K. Leung, "Texture Synthesis by Non-parametric
%  Sampling", ICCV 99.
%
% Notes:
%
%   MAB: SeedImage, and subsequently synthesised_image should have NaNs at 
%   unfilled pixels not zeros. This is now checked for when
%   synthesised_image is created. 06/03/2007
%
% See also:


% Unpack the arguments:
args = u_packargs(varargin, 0, ... % strict mode
          {...  % The mandatory arguments
		  'SeedImage', ...
		  'FilledImage', ...
		  'WindowSize' ...
          'SampleImage',...
          },... % The optional arguments 
          'Debug', 0,...
          'ErrThreshold', 0.01,...
          'Uint8Mode', 1);

      clear varargin;
      
if isempty(args.SampleImage) && ~sum(args.FilledImage(:))
    error('Either the sample image or filled image must be non-empty');
end

if args.Uint8Mode
    % we need SampleImage and SeedImage to be uint8's because we're assuming pixel values in the
    % range 0-255.
    if ~(isa(args.SeedImage, 'uint8') && isempty(args.SampleImage))
        error('SampleImage and SeedImage need to be of type unit8')
    end
    % convert SampleImage to pixels values in the range 0 to 1
    SeedImage = double(args.SeedImage) ./ 255; % SampleImage is now a double
    
    % convert SampleImage to pixels values in the range 0 to 1
    SampleImage = double(args.SampleImage) ./ 255; % SeedImage is now a double     
    
else
    SeedImage = double(args.SeedImage);
    SampleImage = double(args.SampleImage);   
end

% pre-allocation:
if sum(args.FilledImage(:)) == 0 % if there are no filled pixels
    [args.FilledImage, synthesised_image] = seed_the_image(SampleImage, size(SeedImage)); % seed it
else
    synthesised_image = SeedImage;
    synthesised_image(~args.FilledImage) = NaN; %see MAB note in header
end

% make the weighting function -- done here as computing the Gaussian is slow
G = fspecial('gaussian',[args.WindowSize args.WindowSize], args.WindowSize/6.4); % 6.4 as was used by Efros

% Set up visual output for args.Debugging
if args.Debug == 1
    figure;
    subplot(2,3,1);
end

% while there are unfilled pixels in the image
while 1==1  
    % get the list of unfilled pixels
    unfilled_list = get_unfilled_perimeter_pixels(args.FilledImage);
    
    % if there are no more pixels to fill in, then we can stop
    if isempty(unfilled_list)
        break;
    end
    
    % randomly permute the unfilled list
    unfilled_list = unfilled_list(randperm(length(unfilled_list)));
  
    display('Synthesising pixels for next onion skin...');
    for i = 1 :length(unfilled_list)

        % get a window aound the current unfilled_pixel in our (not quite) synthesised_image
        [row, col] = ind2sub(size(synthesised_image), unfilled_list(i));
        sampled_window = sample_window(synthesised_image, args.WindowSize, row, col);
       
        [pixel_list, similarity_image] = find_matches(sampled_window, SampleImage, G, args.ErrThreshold, row, col);
%         [pixel_list, similarity_image] = find_matches2(sampled_window, SampleImage, G, args.ErrThreshold);
        rn = ceil(length(pixel_list) * rand(1));
        % now get the corresponding pixel value from SampleImage
        p = SampleImage(pixel_list(rn));       
        % now we have out p, we can place it in the image
        synthesised_image(row, col) = p;
        args.FilledImage(row, col) = 1;
    end
end

if args.Uint8Mode
    % convert the image to a uint8, in the range 0 to 255
    synthesised_image = uint8(synthesised_image .* 255);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [unfilled_list] = get_unfilled_perimeter_pixels(FilledImage)
%
% Get a list of unfilled pixels neighbouring the pixels already filled.
% This is returned as a list of indices -- use ind2sub to get the image
% locations back.

unfilled_image = bwmorph(FilledImage, 'dilate'); % dilation with SE=ones(3)
unfilled_image = double(unfilled_image) - double(FilledImage);
unfilled_list = find(unfilled_image==1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [FilledImage, synthesised_image] = seed_the_image(SampleImage, sz)
%
% This function seeds the image from the sample image with a 3x3 sample
% taken randomly from the sample image.

% pre-allocate
%synthesised_image = zeros(sz);
synthesised_image = repmat(NaN, sz);
FilledImage = zeros(sz);

% first get the sample 3x3 seed
cont = 1; % dummy variable to replicate a do-while structure
while cont == 1
    row_rand = floor(size(SampleImage,1) * rand) + 1;
    col_rand = floor(size(SampleImage,2) * rand) + 1;
    if row_rand-1 >= 1 && row_rand+1 <= size(SampleImage,1) && ...
            col_rand-1 >= 1 && col_rand+1 <= size(SampleImage,2)
        sample = SampleImage(row_rand-1 : row_rand + 1, col_rand-1 : col_rand+1);
        cont = 0; % OK, we can stop now!
    end
end

% now need to work out where to put it!
% for now, just put it slap-bang in the middle of our image to sample
row_pos = ceil(sz(1) / 2);
col_pos = ceil(sz(2) / 2);

% check that we can put it here and we haven't got a stupidly-sized image
if row_pos-1 >= 1 && row_pos+1 <= sz(1) && ...
        col_pos-1 >= 1 && col_pos+1 <= sz(2)
    synthesised_image(row_pos-1 : row_pos+1 , col_pos-1 : col_pos+1) = sample;
    FilledImage(row_pos-1 : row_pos+1 , col_pos-1 : col_pos+1) = ones(size(sample));
    % this records where we have populated the synthetic image
    disp('Placing the seed in the centre, and this code assumes the seed is 3x3.')
else
    error('We have a problem!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%