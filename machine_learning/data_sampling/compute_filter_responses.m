function [responses] = compute_filter_responses(image_in, decomposition_args)
% Given an image and some arguments, sample feature vectors

f_debug = (nargin == 0 && nargout == 0);
if f_debug, test_script(); return; end

[responses] = func(image_in, decomposition_args);

%% The function
function responses = func(image_in, d_args)

if d_args.normalise
    im_mean = mean(image_in(:));
    im_std = std(image_in(:));
    image_in = (image_in - im_mean) / im_std;
end

responses = [];

%Deal with old cases where the decomp type was just a character string,
%rather than a cell
if ischar(d_args.decomp_type)
    d_args.decomp_type = {d_args.decomp_type};
end

%Now we can switch on the decomp type
for ii = 1:length(d_args.decomp_type)
    switch d_args.decomp_type{ii}

        case {'dt'}
            responses.dt = compute_dual_tree(image_in, d_args.levels(end), d_args.use_nag);
            
        case {'gabor'}
            responses.gabor = compute_gabor_responses(image_in, ...
                                                      d_args.sigma_range,...
                                                      d_args.num_angles);                                                
        case {'gabori'}
            responses.gabor = compute_gabor_responses_d(image_in, ...
                                                      d_args.sigma_range(1),...
                                                      d_args.sigma_range(2),...
                                                      d_args.num_angles);
        case {'g'}
            responses.g = compute_gaussian_responses(image_in, d_args.sigma_range);

        case {'g1d'}
            responses.g1d = compute_gaussian_1st_derivatives(image_in, d_args.sigma_range);
            
        case {'g2d'}
            responses.g2d = compute_gaussian_2nd_derivatives(image_in, d_args.sigma_range);

        case {'g2da'}
            separable_responses = compute_gaussian_2nd_derivatives(image_in, d_args.sigma_range);
            responses.g2d = steer_gaussian_2nd_derivatives(separable_responses, [], d_args.num_angles);
        
        case {'g2di'}
            responses.g2d = ...
                compute_gaussian_2nd_derivatives_d(image_in, ...
                                                   d_args.sigma_range(1), ...
                                                   d_args.sigma_range(2));
                                               
        case {'g2dia'}
            separable_responses = ...
                compute_gaussian_2nd_derivatives_d(image_in, ...
                                                   d_args.sigma_range(1), ...
                                                   d_args.sigma_range(2));
            responses.g2d = steer_gaussian_2nd_derivatives(separable_responses, [], d_args.num_angles);
            
        case {'h2d'}
            responses.h2d = compute_hilbert_2nd_derivatives_sep(image_in, d_args.sigma_range);
            
        case {'h2di'}
            responses.h2d = ...
                compute_hilbert_2nd_derivatives_d(image_in, ...
                                                   d_args.sigma_range(1), ...
                                                   d_args.sigma_range(2));     
        case {'h2da'}
            separable_responses = compute_hilbert_2nd_derivatives_sep(image_in, d_args.sigma_range);
            responses.h2d = steer_hilbert_2nd_derivatives(separable_responses, [], d_args.num_angles);
            
        case {'h2dia'}
            separable_responses = ...
                compute_hilbert_2nd_derivatives_d(image_in, ...
                                                   d_args.sigma_range(1), ...
                                                   d_args.sigma_range(2));
            responses.h2d = steer_hilbert_2nd_derivatives(separable_responses, [], d_args.num_angles);

        case {'haar'}
            responses.haar = compute_haar_responses(image_in, d_args.sigma_range);

        case {'mono'}
            [local_amp local_phase local_ori] = ...
                monogenic(image_in, d_args.num_levels, ...
                                    d_args.min_wavelength, 2, ...
                                    d_args.onf, 1);
            responses.mono = cat(4, local_amp, local_phase, local_ori);
            responses.mono = responses.mono(:,:,2:end,:);
            
        case {'linop'}
            % Do nothing - it's all done in sample_image_features()
            responses.raw = image_in;

        case {'pixel'}
            responses.raw = image_in;
    end
end


%% Old function (for reference)
function features = old_func(channel, rows, cols, decomposition_args)
% Given an image and some arguments, sample feature vectors

switch decomposition_args.decomp_type
    case 'dt' 
        %Compute dual-tree transform of image
        dt = compute_dual_tree(channel, decomposition_args.levels(end), decomposition_args.use_nag);

        %Sample DT coefficients from specified rows and cols according to
        %sampling arguments
        features = sample_dt_data(dt, rows, cols, decomposition_args);
        clear dt;

    case 'dtg2'
        % Sample DT coefficients from specified rows and cols according to
        % sampling arguments
        dt = compute_dual_tree(channel, decomposition_args.levels(end), decomposition_args.use_nag);
        dt_features = sample_dt_data(dt, rows, cols, decomposition_args);
        clear dt;
        
        g2d_responses = compute_gaussian_2nd_derivatives(channel, decomposition_args.sigma_range);
        clear channel;
        
        g2d_features = sample_g2d_data(g2d_responses(:,:,:,1),... 
                                       g2d_responses(:,:,:,2),...
                                       g2d_responses(:,:,:,3), ...
                                       rows, cols, decomposition_args);
        clear g2d_responses;
        
        features = [dt_features g2d_features];
        clear dt_features g2d_features;

    case 'mono'
        [local_amp local_phase local_ori] = ...
            monogenic(channel, decomposition_args.num_levels, decomposition_args.min_wavelength, 2, decomposition_args.onf, 1);
        clear channel;
        
        features = sample_monogenic_data(local_amp, local_phase, local_ori, rows, cols, decomposition_args);
        clear local_amp local_phase;

    case 'g1d'
        g1d_responses = compute_gaussian_1st_derivatives(channel, decomposition_args.sigma_range);
        clear channel;
        
        g1d_features = sample_g1d_data(g1d_responses(:,:,:,1),... 
                                       g1d_responses(:,:,:,2),...
                                       rows, cols, decomposition_args);
        clear g1d_responses;

        features = g1d_features;
        clear g1d_features;

    case 'g2d'
        g2d_responses = compute_gaussian_2nd_derivatives(channel, decomposition_args.sigma_range);
        clear channel;
        
        g2d_features = sample_g2d_data(g2d_responses(:,:,:,1),... 
                                       g2d_responses(:,:,:,2),...
                                       g2d_responses(:,:,:,3), ...
                                       rows, cols, decomposition_args);
        clear g2d_responses;
        
        features = g2d_features;
        clear g2d_features;

    case 'g2di'       
        g2d_responses = compute_gaussian_2nd_derivatives_d(channel, decomposition_args.sigma_range(1), decomposition_args.sigma_range(2));
        clear channel;
        features = sample_g2d_data_d(g2d_responses, rows, cols, decomposition_args);
        clear g2d_responses;

    case 'g12d'
        g1d_responses = compute_gaussian_1st_derivatives(channel, decomposition_args.sigma_range);
        g2d_responses = compute_gaussian_2nd_derivatives(channel, decomposition_args.sigma_range);
        clear channel;

        g1d_features = sample_g1d_data(g1d_responses(:,:,:,1),... 
                                       g1d_responses(:,:,:,2),...
                                       rows, cols, decomposition_args);
        clear g1d_responses;

        g2d_features = sample_g2d_data(g2d_responses(:,:,:,1),... 
                                       g2d_responses(:,:,:,2),...
                                       g2d_responses(:,:,:,3), ...
                                       rows, cols, decomposition_args);
        clear g2d_responses;
        
        features = [g1d_features g2d_features];
        clear g2d_features g1d_features;

    case 'clover'
        g2d_responses = compute_clover_responses(channel, decomposition_args.sigma_range);
        clear channel;
        
        g2d_features = sample_g2d_data(g2d_responses(:,:,:,1),... 
                                       g2d_responses(:,:,:,2),...
                                       g2d_responses(:,:,:,3), ...
                                       rows, cols, decomposition_args);
        clear g2d_responses;
        
        features = g2d_features;
        clear g2d_features;

    case 'haar'
        haar_responses = compute_haar_responses(channel, decomposition_args.sigma_range);
        clear channel;
        
        haar_features = sample_haar_data(haar_responses(:,:,:,1),... 
                                         haar_responses(:,:,:,2),...
                                         rows, cols, decomposition_args);
        clear haar_responses; 
        
        features = haar_features;
        clear haar_features;

    case 'pixel'
        %Sample pixels from image
        features = sample_pixel_data(channel, rows, cols, decomposition_args);

    case 'linop'
        features = sample_linop_data(channel, rows, cols, decomposition_args);
        clear channel;

    otherwise
        if strcmp(get_username,'ptresadern')
            error(['Unknown decomposition: ',decomposition_args.decomp_type]);
        end
end


%% Test script
function test_script()
clc; 

rand('twister',3141592);
image_in = round(rand(256,256)*255);

inds = ceil(rand(100,1)*numel(image_in));
[rows,cols] = ind2sub(size(image_in), inds);
rows = min(max(rows, 2), 254);
cols = min(max(cols, 2), 254);

args = default_args();
args.num_levels = 4;
args.subtract_mean = true;

for decomp_type = {'pixel', ...
                   'linop', 'dt', 'dtg2', 'mono', ...
                   'g1d', 'g2d', 'clover', 'haar', 'g2di', 'g12d' }
    args.decomp_type = decomp_type{1};
    decomposition_args = get_decomposition_args_from(args);

    fprintf('Testing %s...', decomp_type{1});
    
    warning('off', 'ASYM:unexpectedArgument');
    old_features = old_func(image_in, rows, cols, decomposition_args);
    responses = func(image_in, decomposition_args);
    new_features = sample_image_features(responses, rows, cols, ...
                                         decomposition_args);
    warning('on', 'ASYM:unexpectedArgument');
   
    if all(old_features(:) == new_features(:))
        fprintf('success\n');
    else
        fprintf('fail\n');
    end
end
