function compare_g2d_dt()

% clc; 
clear;

warning('off', 'ASYM:unexpectedArgument');

rand('state',3141592);
randn('state',3141592);

line_width = 6;

bgc = line_width*6;
bgsz = 2*bgc;

do_g2d = true;
do_dt = false;

do_linear = true;
do_tree = false;
do_ann = false;

normalize_inputs = false;

if do_dt,   dt_levels = 1:3;
else        dt_levels = [];
end
n_dt_levels = length(dt_levels);

n_sigma = n_dt_levels*12 / 3;
% sigma_range = 2.^linspace(0, log2(4), n_sigma);
sigma_range = [2.5]; n_sigma = length(sigma_range); n_params = n_sigma*3;

contrast = 8;
noise_sigma = contrast * 0.1;

xo_range = 0;
yo_range = xo_range;

% generate 1000 example images of lines with a fixed width and no noise
args = u_packargs([], ... % the user's input
    '0', ... % non-strict mode
    'bg_size',				unixenv('BG_SIZE',[bgsz bgsz]), ...
    'bg_type',				unixenv('BG_TYPE','flat'), ...
    'orientation_range',    unixenv('ORIENTATION_RANGE',[0 180]), ...
    'width_range', 			unixenv('WIDTH_RANGE',line_width*[0.25 1]), ...
    'contrast_range',       unixenv('CONTRAST_RANGE',contrast*[1 1]), ...
    'squash_range',         unixenv('SQUASH_RANGE',[0 0]), ...
    'line_type',			unixenv('BAR_TYPE','sin'), ...
    'win_size',				unixenv('WIN_SIZE',1), ...
    'rgb_channel',          unixenv('RGB_CHANNEL','mono'), ...
    'use_nag',				unixenv('USE_NAG',false),...
    'normalise', 			unixenv('NORMALISE',0), ...
    ...
    'levels',               unixenv('LEVELS',dt_levels), ...
    'do_max',				unixenv('DO_MAX',false), ...
    'rotate',				unixenv('ROTATE',false), ...
    'decomp_type', 			unixenv('DECOMP_TYPE','g2d'), ...
        ... % DTCWT parameters
        'feature_shape', 		unixenv('FEATURE_SHAPE','rect'), ...
        'feature_type',			unixenv('FEATURE_TYPE','conj'), ...
        ... % Gaussian derivative parameters
        'sigma_range', 			unixenv('SIGMA_RANGE',sigma_range), ...
    'prediction_type',		unixenv('PREDICTION_TYPE','rf_regression'), ...
        ... % Tree/Forest parameters
        'n_trees',				unixenv('NUM_TREES',200), ...
        'split_criterion_c',    unixenv('SPLIT_CRITERION_C','gdi'),...
        'var_criterion_c',		unixenv('VAR_CRITERION_C','mabs'),...
        'split_criterion_r',    unixenv('SPLIT_CRITERION_R','dabs'),...
        'var_criterion_r',		unixenv('VAR_CRITERION_R','mabs'),...
        'split_min',			unixenv('SPLIT_MIN',10), ...
        'end_cut_min',			unixenv('END_CUT_MIN',1), ...
        'do_ubound',			unixenv('DO_UBOUND',1), ...
        'do_circular',			unixenv('DO_CIRCULAR',true), ...
        'w_prior',				unixenv('W_PRIOR',0), ...
        'impure_thresh',		unixenv('IMPURE_THRESH',1e-4), ...
        'minimise_size',		unixenv('MINIMIZE_TREE',0), ...
        'd',                    unixenv('d',[]), ...
        'quiet',                unixenv('QUIET', 0), ...
        ... % Neural network parameters
        'n_iterations',         unixenv('ANN_N_ITERATIONS', 100), ...
        'verbose',              unixenv('ANN_VERBOSE', false), ...
    'dummy',                [] ...
);

n_params = max(n_sigma*3, n_dt_levels*12) * args.win_size;
n_params = 3;
n_images = n_params * 25; display(n_images);
% n_images = 5000;

args.num_levels = length(args.levels);
args.random_m = ceil(n_params/2);

if do_g2d
    args.decomp_type = 'g2d';
    g2d_args = get_decomposition_args_from(args);
    g2d_samples = get_samples_per_channel(g2d_args);
end

if do_dt
    args.decomp_type = 'dt';
    dt_args = get_decomposition_args_from(args);
    dt_samples = get_samples_per_channel(dt_args);
end

% Generate a line aligned with the filters to compute maximum responses
if do_g2d
    tmpargs = args; 
    tmpargs.orientation_range = [0 180];
    tmpargs.width_range = line_width*[1 1];
    [im, lbl, lblc, lblo, parameters] = create_line_image(tmpargs);
    g2d_responses = compute_filter_responses(im, g2d_args);
    g2d_maxresp = sample_image_features(g2d_responses, bgc, bgc, g2d_args);

    if (args.win_size == 1)
        r = reshape(g2d_maxresp,[n_sigma,3]);
        g2d_maxresp = r(:,3);
        [ignore,maxind] = min(g2d_maxresp);
        fprintf('Optimal sigma = %.2f\n', sigma_range(maxind));
    end
    
%     filts = gabor_filters(90, 4.0);
%     hw = (size(filts,1)-1) / 2;
%     xoff = 1;
%     im = im(bgc+xoff-hw:bgc+xoff+hw, bgc-hw:bgc+hw);
%     r = zeros(size(filts,3),2);
%     for t = 1:size(filts,3)
%         r(t,1) = conv2(real(filts(:,:,t)), im, 'valid');
%         r(t,2) = conv2(imag(filts(:,:,t)), im, 'valid');
%     end
%     parameters.orientation
%     figure(5); clf; hold on;
%     x = linspace(0,360,size(filts,3)*2+1);
%     x = x(1:end-1);
%         plot(x, [r(:,1); r(:,1)], 'b-');
%         plot(x, [r(:,2); -r(:,2)], 'r-');
%     maxr = max(abs(r(:,1))); minr = min(abs(r(:,1)));
%     x = linspace(-90,90,size(r,1));
%     sig = 24;
%     y = minr + (maxr-minr) * exp(-0.5 * x.*x / (sig*sig));
%         plot(x+parameters.orientation+90, y, 'g-');
end

if ~exist('g2d_inputs','var')
    if do_g2d
        g2d_inputs = zeros(n_images, g2d_samples);
    end
    if do_dt
        dt_inputs = zeros(n_images, dt_samples);
    end

    targets = zeros(n_images, 1);
    widths = zeros(n_images, 1);

    fprintf('Generating training data...');
    for i = 1:n_images
        [im, lbl, lblc, lblo, parameters] = create_line_image(args);
        im = im + randn(size(im))*noise_sigma;

        theta = parameters.orientation;
        c2t = cosd(2*theta); s2t = sind(2*theta);
        targets(i) = complex(c2t, s2t);
        widths(i) = parameters.width;

        xo = round((2*rand-1)*xo_range);
        yo = round((2*rand-1)*yo_range);

        % for each image, compute the g2d and dt responses for the centre pixel
        if do_g2d
            g2d_responses = compute_filter_responses(im, g2d_args);
            g2d_inputs(i, 1:g2d_samples) = ...
                sample_image_features(g2d_responses, bgc+yo, bgc+xo, g2d_args);
        end

        if do_dt
            dt_responses = compute_filter_responses(im, dt_args);
            dt_inputs(i, 1:dt_samples) = ...
                sample_image_features(dt_responses, bgc+yo, bgc+xo, dt_args);
        end
    end
    fprintf('done\n');
end

if normalize_inputs
    g2d_inputs = normalize_g2d(g2d_inputs, n_images, n_sigma);
end

theta_gt = angle(targets);

fprintf('Building predictors...');
if do_g2d
    if (args.win_size == 1)
        cols = 1:n_sigma;       Ixy = g2d_inputs(:,cols);
        cols = cols + n_sigma;  Ixx = g2d_inputs(:,cols);
        cols = cols + n_sigma;  Iyy = g2d_inputs(:,cols);
        A = (Ixx+Iyy) / 2;
        B = sqrt(((Ixx-Iyy) / 2).^2 + Ixy.^2);
        g2d_resp_power = A-B;

        figure(1); clf; hold on;
            plot(sigma_range, g2d_resp_power, '.', ...
                 'color', 0.75*ones(1,3));
            plot(sigma_range, g2d_maxresp, '.', ...
                 'color', 'b');
            plot(sigma_range(maxind), g2d_maxresp(maxind), '.', ...
                 'color', 'r');
    end
    
    X = g2d_inputs; y = targets;
    
    if do_linear
        % fit a linear model to each of the datasets
        predictor = pt_lin_reg_train(X, real(y));
        predictor_imag = pt_lin_reg_train(X, imag(y));
        predictor.beta = complex(predictor.beta, predictor_imag.beta);
        g2d_pred = predictor;

        if (args.win_size == 1)
            b = reshape(g2d_pred.beta(2:end),[n_sigma,3]);
%             display(b);
            sigma = sigma_range(n_sigma);
            [g,dg,ddg] = gaussian_filters_1d(sigma);
            hw = (length(g)-1)/2;
            Freal = zeros(2*hw+1); Fimag = zeros(2*hw+1);
            for i = 1:n_sigma
                sigma = sigma_range(i);
                [g,dg,ddg] = gaussian_filters_1d(sigma, hw);
                Gxy = -dg'*dg;
                Gxx = g'*ddg;
                Gyy = ddg'*g;
                Freal = Freal + ...
                        real(b(i,1))*Gxy + ...
                        real(b(i,2))*Gxx + ...
                        real(b(i,3))*Gyy;
                Fimag = Fimag + ...
                        imag(b(i,1))*Gxy + ...
                        imag(b(i,2))*Gxx + ...
                        imag(b(i,3))*Gyy;
            end
            figure(10); clf; colormap(gray(256));
                imagesc([Freal Fimag]); axis('image','off');
%             display([min(Freal(:)) max(Freal(:)); min(Fimag(:)) max(Fimag(:))]);
        
            g = reshape(g2d_inputs,[n_images,n_sigma,3]);
            s = mean(sqrt(sum(g.^2,3)));

            b = g2d_pred.beta;
            b(1,:) = [];
            b = reshape(b,[n_sigma,3])';
            d = sqrt(diag(b'*b));

            % Show distribution of weights wrt sigma
            figure(2); clf; hold on;
                plot(sigma_range, d/max(d), 'b.-');
                plot(sigma_range, s, 'r.-');

            % ideal values
            g2d_analytic = g2d_pred;
            g2d_analytic.beta = zeros(n_sigma,3);
            g2d_analytic.beta(maxind,:) = ...
                [complex(0,2) complex(1,0) complex(-1,0)];
            g2d_analytic.beta = [0; g2d_analytic.beta(:)];
        end
    end
    
    if do_tree
        g2d_tree = cell(1,args.n_trees);
        for i = 1:args.n_trees
            rows = randperm(n_images);
            rows = rows(1:ceil(n_images/2));
            g2d_tree{i} = tree_reg_train(X(rows,:), y(rows), args);
        end
    end
    
    if do_ann
        g2d_ann = pt_ann_reg_train(X, [real(y) imag(y)], args);
    end
end

if do_dt
    X = dt_inputs; y = targets;
    
    if do_linear
        predictor = pt_lin_reg_train(X, real(y));
        predictor_imag = pt_lin_reg_train(X, imag(y));
        predictor.beta = complex(predictor.beta, predictor_imag.beta);
        dt_pred = predictor;
    end
    
    if do_tree
        dt_tree = cell(1,args.n_trees);
        for i = 1:args.n_trees
            rows = randperm(n_images);
            rows = rows(1:ceil(n_images/2));
            dt_tree{i} = tree_reg_train(X(rows,:), y(rows), args);
        end
    end

    if do_ann
        dt_ann = pt_ann_reg_train(X, [real(y) imag(y)], args);
    end
end
fprintf('done\n');

rand('state',27183);
randn('state',27183);

% generate another 1000 images and compute prediction errors
fprintf('Generating test data...');
for i = 1:n_images
    [im, lbl, lblc, lblo, parameters] = create_line_image(args);
    im = im + randn(size(im))*noise_sigma;
    
    theta = lblo(bgc,bgc);
    c2t = cosd(2*theta); s2t = sind(2*theta);
    targets(i) = complex(c2t, s2t);

    xo = round((2*rand-1)*xo_range);
    yo = round((2*rand-1)*yo_range);

    if do_g2d
        g2d_responses = compute_filter_responses(im, g2d_args);
        g2d_inputs(i, 1:g2d_samples) = ...
            sample_image_features(g2d_responses, bgc+yo, bgc+xo, g2d_args);
    end
    
    if do_dt
        dt_responses = compute_filter_responses(im, dt_args);
        dt_inputs(i, 1:dt_samples) = ...
            sample_image_features(dt_responses, bgc+yo, bgc+xo, dt_args);
    end
end
fprintf('done\n');

if normalize_inputs
    g2d_inputs = normalize_g2d(g2d_inputs, n_images, n_sigma);
end

if do_g2d
    X = g2d_inputs; y = targets;
    
    if do_linear
        if exist('g2d_analytic','var')
            y_fit = linear_regressor_predict(g2d_analytic, X);
            o_errors = ori_error(y_fit, targets) * 180/pi;
            fprintf('%20s: ', 'g2d analytic');
            disp([mean(abs(o_errors)) std(o_errors)]);
        end
        
        y_fit = linear_regressor_predict(g2d_pred, X);
        o_errors = ori_error(y_fit, targets) * 180/pi;
        fprintf('%20s: ', 'g2d linear');
        disp([mean(abs(o_errors)) std(o_errors)]);
    end
    
    if do_tree
        y_fit = tree_predict(g2d_tree{1}, X);
        for i = 2:args.n_trees
            y_fit = y_fit + tree_predict(g2d_tree{i}, X);
        end
        y_fit = y_fit / args.n_trees;
        o_errors = ori_error(y_fit, targets) * 180/pi;
        fprintf('%20s: ', 'g2d forest');
        disp([mean(abs(o_errors)) std(o_errors)]);
    end
    
    if do_ann
        y_fit = ann_regressor_predict(g2d_ann, X);
        y_fit = complex(y_fit(:,1), y_fit(:,2));
        o_errors = ori_error(y_fit, targets) * 180/pi;
        fprintf('%20s: ', 'g2d ann');
        disp([mean(abs(o_errors)) std(o_errors)]);
    end  
end

if do_dt
    X = dt_inputs; y = targets;
    
    if do_linear
        y_fit = linear_regressor_predict(dt_pred, X);
        o_errors = ori_error(y_fit, targets) * 180/pi;
        fprintf('%20s: ', 'dt linear');
        disp([mean(abs(o_errors)) std(o_errors)]);
    end

    if do_tree
        y_fit = tree_predict(dt_tree{1}, X);
        for i = 2:args.n_trees
            y_fit = y_fit + tree_predict(dt_tree{i}, X);
        end
        y_fit = y_fit / args.n_trees;
        o_errors = ori_error(y_fit, targets) * 180/pi;
        fprintf('%20s: ', 'dt forest');
        disp([mean(abs(o_errors)) std(o_errors)]);
    end
    
    if do_ann
        y_fit = ann_regressor_predict(dt_ann, X);
        y_fit = complex(y_fit(:,1), y_fit(:,2));
        o_errors = ori_error(y_fit, targets) * 180/pi;
        fprintf('%20s: ', 'dt ann');
        disp([mean(abs(o_errors)) std(o_errors)]);
    end
    
    if (args.win_size == 1)
        % Show distribution of weights wrt levels
        m = reshape(mean(dt_inputs),[6,n_dt_levels*2])';
        display(m);

        b = dt_pred.beta;
        b(1,:) = [];
        d = abs(b);
        d = reshape(d,[args.num_levels*6,2])';
        d = reshape(d,[12,args.num_levels]);
        figure(3); clf;
            plot(d(1:2:end,:));
    end
end

warning('on', 'ASYM:unexpectedArgument');


function [g2d_out] = normalize_g2d(g2d_in, n_images, n_sigma)
g2d_in = reshape(g2d_in, [n_images,n_sigma,3]);
g2d_out = zeros(size(g2d_in));
g2d_out(:,1:end-1,:) = g2d_in(:,2:end,:) - g2d_in(:,1:end-1,:);
g2d_out(:,end,:) = g2d_in(:,1,:);
g2d_out = g2d_out(:,:);

