function [fg_idx, bg_idx] = sample_from_maps(fg_map, bg_map, ...
                                             samples_per_image, sampling_args)
% Given a foreground and background probability map, select pixel positions for
% foreground and background samples

% Map can be either binary mask (sample uniformly from the set, without 
% replacement) or a real-valued probability map (use weighted sampling).

% Compute proportions of foreground and background pixels
% Ratio of fg:bg is proportional to 1:args.bg_ratio
fg_prop =                      1 / (1 + sampling_args.bg_ratio);
bg_prop = sampling_args.bg_ratio / (1 + sampling_args.bg_ratio);

% Q: What if samples_per_image > #pixels in foreground/background?

% Select foreground sample positions (all methods)
switch sampling_args.output_type
    case {'dummy'}
        % Balance sampling so that samples are roughly equally distributed
        % across vessel widths

    otherwise
        % Uniform sampling
      
        % Get random sample of vessel pixels
        n_samples = floor(samples_per_image * fg_prop);
        fg_idx = sample_from_map(fg_map, n_samples, sampling_args.replace_sample);
end            

% Choose background pixels (detection methods only)
if bg_prop > 0

    % Get random sample of background pixels
    n_samples = ceil(samples_per_image * bg_prop);
    bg_idx = sample_from_map(bg_map, n_samples, sampling_args.replace_sample);
        
else
    bg_idx = [];
end


function idx = sample_from_map(map, n_samples, replace_sample)
% Sample points from an image map

if islogical(map)
    % Binary map - sample uniformly
    idx = find(map);
    
    %with replacement
    if replace_sample
        idx = ceil(numel(idx)*rand(n_samples,1));
    else
    %without replacement
        if numel(idx) > n_samples
            shuffled = randperm(numel(idx));
            idx = idx(shuffled(1:n_samples));
        else
            display('Warning, not enough samples in image');
        end
    end
else
    % Real-valued probabilities - sample in proportion to weights
    idx = sample_from_probs(map, n_samples);
    
%     %DEGUG CODE
%     [y x] = ind2sub(size(map), idx);
%     figure; imgray(map); plot(x, y, 'rx');
    
end


