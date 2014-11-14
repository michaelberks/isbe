function [a_shapes, a_scales, mean_target, a_rots, a_trans a_origins]...
    = align_shapes(u_shapes, varargin)
% ALIGN_SHAPES  
%    [a_shapes] = align_shapes(u_shapes, target_area, ...)
%
%    inputs:
%       u_shapes  - N x 2*dim matrix of unaligned shapes
%      optional:
%       seed      - {1} Index of initial shape to align to
%       target_area      - target_area of mean (in pixels) to align to
%       setOrigin - {0}/1, flag to determine whether origin is optimised
%
%    outputs:
%       a_shapes  - N x 2*dim matrix of aligned shapes
%
%    notes:
%
% Copyright: (C) 2006-2008 University of Manchester
% Author: Michael Berks



args = u_packargs(varargin, 0,...
        'seed', 1,...
        'shiftOrigin', 1,...
        'target_area', 50000,...
        'length', 0,...
        'tol', 1e-4);
    clear varargin
target_area = args.target_area;
target_length = args.length;
seed = args.seed;

%
% first align vectors using procrustes analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

target_shape = reshape(u_shapes(seed, :), [], 2);
%[a_shapes, a_scales, sq_diff] =...
[a_shapes, sq_diff] =...
    pro_group(target_shape, u_shapes([1:seed-1,seed+1:end],:), 0);

a_shapes = [a_shapes(1:seed-1,:); target_shape(:)';...
            a_shapes(seed:end,:)];
% a_scales = [a_scales(1:seed-1); 1; a_scales(seed:end)];

clear target_shape;


mean_sq_diff = mean(sq_diff);
display(['mean change = ', num2str(mean_sq_diff)]);

go_on = 1; jj = 1;
while go_on
    
    mean_shape = mean(a_shapes);
    mean_shape = reshape(mean_shape, [], 2);
    mean_shape = mean_shape - repmat(mean(mean_shape), size(mean_shape,1), 1);
    
    if target_area
        current_area = polyarea(mean_shape(:,1), mean_shape(:,2));
        scale_factor = sqrt(target_area / current_area);
    else
        current_length = sum(sqrt(sum(diff(mean_shape).^2,2)));
        scale_factor = target_length / current_length;
    end
    mean_shape = mean_shape*scale_factor;

    th = atan2(mean_shape(1,1), mean_shape(1,2));
    mean_shape = ([cos(th) -sin(th); sin(th) cos(th)] * mean_shape')';
    
    [a_shapes, sq_diff] =...
        pro_group(mean_shape, a_shapes, args.shiftOrigin);
%     a_scales = a_scales.* it_scales;
    
    go_on = abs(mean_sq_diff - mean(sq_diff)) > args.tol;
    mean_sq_diff = mean(sq_diff);
    display(['mean change = ', num2str(mean_sq_diff)]);
    
    jj = jj+1;
end
mean_target = mean_shape;
display(['stopped after ', num2str(jj-1), ' iterations']);

% Perform final alignment on group from original unaligned shapes to final
% target shape to obtain pose parameters. Would be much nicer and more
% efficient to track these as we go, but haven't implemented this
[a_shapes dummy a_scales a_rots a_trans a_origins] =...
                    pro_group(mean_shape, u_shapes, args.shiftOrigin);

end

    function [aligned_shapes diffs scales rots trans origins] =...
                    pro_group(target_shape, shapes, shiftOrigin)
    
        [N dim] = size(shapes); dim = dim/2; % N is number of shape vectors
        %dim is number of points in each vector
        aligned_shapes = zeros(N, 2*dim);
        scales = zeros(N, 1);
        rots = zeros(2, 2, N);
        trans = zeros(N, 2);
        diffs = zeros(N, 1);
        origins = ones(N,1);

        for ii=1:N
            source_shape = reshape(shapes(ii, :), [], 2);

            if shiftOrigin
                % cycle around the shape, shifting the origin each time
                % we choose the one with the minimal SSD to the target
                ssd = zeros(dim, 1); 
                for origin = 1 : dim
                    shift_shape = circshift(source_shape, origin - 1); 
                    ssd (origin) = mb_procrustes(target_shape, shift_shape, []);
                    %ssd (origin) = procrustes(target_shape, shift_shape, 'Reflection', 'best');
                end
                [dummy, origin] = min (ssd);
                origins(ii) = origin;
                % apply the optimal origin shift to source shape...
                source_shape = circshift(source_shape, origin - 1);
            end

            [dd Z t] = mb_procrustes(target_shape, source_shape, []);
            %[dd Z t] = procrustes(target_shape, source_shape, 'Reflection', 'best');
            aligned_shapes(ii, :) = [Z(:,1)', Z(:,2)'];
            scales(ii) = t.b;
            rots(:,:,ii) = t.T;
            trans(ii,:) = t.c(1,:) / t.b;
            diffs(ii) = dd;
            %%%
            % pose relationship between source and target is
            % tgt = scale*(src*rot - trans) 
        end
    end
        
