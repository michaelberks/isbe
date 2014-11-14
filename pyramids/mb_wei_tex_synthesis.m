function [synthesised_image ] = mb_wei_tex_synthesis(varargin)
%
% MB_WEI_TEX_SYNTHESIS Synthesise an image using an ergodic GMM.
%
% Mandatory Arguments:
% 'TargetImage'
%  - Hello
%
% 'FilledImage'
% -
%
% 'SampleImage1'
%   -
%
% Optional Arguments:
%
% 'SampleImage2'
%   - 
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
% References:
%
%
% See also:
%

% Unpack the arguments:
args = u_packargs(varargin, 0, ... % strict mode
		  {...  % The mandatory arguments
          'TargetImage1',...
          'TargetImage2'...
          },... % The optional arguments
          'PixelOrder', [],...
          'FilledImage', [],...
          'SampleImage1', [], ...
          'SampleImage2', [], ...
          'SampleData', [], ...
          'BothLevels', 0,...
          'SaveFrequency', 1000,...
          'SaveFile', [], ...
          'WindowSize1', 5,...
          'WindowSize2', 5,...
          'Quiet', true,...
          'K_method', 'standard');

synthesised_image = args.TargetImage2;
args = rmfield(args, 'TargetImage2');
%synthesised_image(~args.FilledImage) = NaN;  

if ~isempty(args.SampleData)
    %Don't need to construct k-nearest data structure
    knearest_data = args.SampleData;
    args = rmfield(args, 'SampleData');
    
elseif ~isempty(args.SampleImage1) && ~isempty(args.SampleImage2)
    
    %Construct k-nearest data structure from sample images
    
    k_args.Image1 = args.SampleImage1;
    k_args.Image2 = args.SampleImage2;
    k_args.WindowSize1 = args.WindowSize1;
    k_args.WindowSize2 = args.WindowSize2;
    k_args.Method = args.K_method;
    k_args.BothLevels = args.BothLevels;

    knearest_data = mb_knearest_build(k_args);
    
    clear k_args;
 
else
    error('Must supply either sample data or images');
end

if isempty(args.PixelOrder)
    if ~isempty(args.FilledImage)
        args.PixelOrder = mb_compute_pixel_order(args.FilledImage);
    else
        error('Must supply either pixel order or filled image');
    end
end

[rows2, cols2] = ind2sub(size(synthesised_image), args.PixelOrder);
rows1 = ceil(rows2/2);
cols1 = ceil(cols2/2);

%We do some manual progress feedback to the user, and set it up
%here. First, define a set of percentages to DISP at.
progress_report_percentages = 0.1:0.1:1;
pixels_to_fill = length(args.PixelOrder);
pixels_filled = 0;
tic; % Start a record of when we started.

% Here is the main algorithm.

for ii = 1:pixels_to_fill
    
%     % choose one of the unfilled_list at random
%     unfilled_pixel_chosen = args.PixelOrder(ii);
%     
%     %%%%%%
%     % get a window aound the current unfilled_pixel in our (not quite) synthesised_image
%     [row2, col2] = ind2sub(size(synthesised_image), unfilled_pixel_chosen);
%     row1 = ceil(row2/2);
%     col1 = ceil(col2/2);
    
    row1 = rows1(ii); col1 = cols1(ii);
    row2 = rows2(ii); col2 = cols2(ii);
    
    sampled_window = reshape(sample_window(args.TargetImage1, ...
					  args.WindowSize1, row1, col1), 1, []);
    sampled_window = (sampled_window - knearest_data.Mean1) / knearest_data.Std1;
    
    if args.BothLevels
        sampled_window2 = reshape(sample_window(synthesised_image, ...
					  args.WindowSize2, row2, col2), 1, []);
        sampled_window2 = (sampled_window2 - knearest_data.Mean2) / knearest_data.Std2;          
        sampled_window = [sampled_window, sampled_window2]; %#ok
    end
    
    % Sample from the sample image(s):
    [values] = mb_knearest_search('Data', knearest_data,...
                                        'TestVector', sampled_window,...
                                        'Method', args.K_method);
    
    if ~args.BothLevels
        %need to choose appropriate value dependent on position relative to
        %lower subscript
        
        values = values(4 - 2*rem(row2, 2) - rem(col2, 2));
    end

    synthesised_image(row2, col2) = values;
    
    % Update the record of where we've been.
    pixels_filled = pixels_filled +  1;
    
    if ~args.Quiet
        % Report progress
        progress_so_far = pixels_filled / pixels_to_fill;
        if progress_so_far > progress_report_percentages(1)

          disp(['Progress: ' num2str(floor(100 * progress_so_far)) ' percent' ...
                ' complete.']);
          disp(['--Time: ' datestr(now)]);
          % Remove the first report point from the progress_report_percentages
          progress_report_percentages(1) = [];
        end

        % see if we need to save
        if ~(mod(pixels_filled, args.SaveFrequency)) && ~isempty(args.SaveFile)
            save(args.SaveFile, 'synthesised_image');
        end
    end
end

% save final image
if ~isempty(args.SaveFile)
    save(args.SaveFile, 'synthesised_image');
end

% Report the generation of the synthetic image
disp('Progress: 100 percent complete.');
disp('A synthetic ergodic texture has been synthesised!');

%Main function end