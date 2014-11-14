function [fg_idx, bg_idx] = sample_vessel_labels(fg_mask, bg_mask, ...
                                                 num_samples_image, args)
% Given a foreground and background mask, select pixel positions for
% foreground and background samples

% Select foreground sample positions (all methods)
switch args.prediction_type
    case {'dummy'}
        % Balance sampling so that samples are roughly equally distributed
        % across vessel widths

    otherwise
        % Uniform sampling
      
        % Get random sample of vessel pixels
        fg_idx = find(fg_mask);
        shuffled = randperm(numel(fg_idx));
        num_fg = floor(num_samples_image / (1 + args.bg_ratio));
        fg_idx = fg_idx(shuffled(1:num_fg));
end            

% Choose background pixels (detection methods only)
switch args.prediction_type
    case {'detection', 'centre_detection', 'foveal_detection', ...
          'logistic_classification', ...
          'junction_detection', 'junction_centre_detection'}
        % Uniform sampling

        % Get random sample of background pixels
        bg_idx = find(bg_mask);
        shuffled = randperm(numel(bg_idx));
        num_bg = ceil(num_samples_image * args.bg_ratio / (1 + args.bg_ratio));
        bg_idx = bg_idx(shuffled(1:num_bg));
        
    otherwise
        bg_idx = [];
end

            
