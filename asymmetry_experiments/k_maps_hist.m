function [var_out] = k_maps_hist(ori_type, line_type, varargin)
%K_MAPS_HIST *Insert a one line summary here*
%   [] = k_maps_hist()
%
% Inputs:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 13-Jun-2011
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    'num_bins', 60,...
    'do_vm', 0,...
    'thresh', 0,...
    'use_lines', 1,...
    'data_type', '\2004_screening_processed\temp_roi\',...
    'ori_stem', 'C:\isbe\asymmetry_project\data\orientation_maps\',...
    'line_stem', 'C:\isbe\asymmetry_project\data\line_maps\',...
    'mam_names', 'C:\isbe\asymmetry_project\data\mam_names\2004_screening_abnormals.mat',...
    'mask_dir', 'C:\isbe\asymmetry_project\data\masks\2004_screening_processed\temp_roi\',...
    'save_path', []);
clear varargin;

%load names of mass regions
abnormal_names = u_load(args.mam_names);

%Count how many regions
num_roi = length(abnormal_names);

%set num_bins
num_bins = args.num_bins;

%pre-allocate outputs
var_out.n_hist_all = zeros(num_roi, num_bins, 2);
var_out.n_hist_line = zeros(num_roi, num_bins, 2);
var_out.N_hist_all = zeros(num_roi, num_bins, 2);
var_out.N_hist_line = zeros(num_roi, num_bins, 2);

%set histograms bins
bins = linspace(-pi/2, pi/2, num_bins);

px_per_mm = 100/9;
r_max = 400;
r_min = px_per_mm*5;
R = px_per_mm*4;

xx = repmat(-400:399, 800, 1);
yy = repmat((399:-1:-400)', 1, 800);

theta = mod(atan2(yy, xx), pi);
rij = sqrt(xx.^2 + yy.^2);
phi = abs(asin(R./rij));
pij = 2*R./(pi*rij);

pvm = pi*R ./ rij;

r_mask = (rij < r_max) & (rij > r_min);

var_out.N_all = sum(r_mask(:));
var_out.N_line = zeros(num_roi, 2);
var_out.p_all = sum(pij(r_mask)) / var_out.N_all;
var_out.p_line = zeros(num_roi, 2);
var_out.pq_line = zeros(num_roi, 2);
var_out.q_line = zeros(num_roi, 2);
var_out.p_line_vm = zeros(num_roi, 2);

for ii = 1:num_roi
    display(['processing region ', num2str(ii)]);
    
    mam_name = cell(2,1);
    mam_name{1} = abnormal_names{ii};
    if abnormal_names{ii}(4) == 'R'
        mam_name{2} = [mam_name{1}(1:3) 'L' mam_name{1}(5:6)];
    else
        mam_name{2} = [mam_name{1}(1:3) 'R' mam_name{1}(5:6)];
    end

    for aa = 1:2

        %Load in regions associated with mass
        mam_mask = u_load([args.mask_dir mam_name{aa} '_roi.mat']);
        
        ori_map = u_load([args.ori_stem ori_type args.data_type ...
            mam_name{aa} '_roi.mat']);
        
        if ~isreal(ori_map)
            ori_pm = ori_map;
            ori_map = angle(ori_pm);
        else
            ori_pm = exp(i*ori_map);
        end
        
        if isempty(line_type)
            line_map = abs(ori_pm);     
        else
            line_map = u_load([args.line_stem line_type args.data_type ...
                mam_name{aa} '_roi.mat']);
        end
        
        if ~isreal(line_map)
            line_map = abs(line_map);
        end
                
        %Work out which pixels point towards the centre
        theta_diff = abs(theta - ori_map);
        theta_mask = ((theta_diff < phi) | (abs(theta_diff-pi) < phi)) & r_mask & mam_mask;

        %Workout which pixels are lines
        line_mask = (line_map > args.thresh) & r_mask & mam_mask;
        
        %Take count of line pixels
        var_out.N_line(ii,aa) = sum(line_mask(:));
        var_out.p_line(ii,aa) = sum(pij(line_mask)) / var_out.N_line(ii,aa);
        var_out.q_line(ii,aa) = sum(pij(~line_mask & r_mask));
        var_out.pq_line(ii,aa) = sum(pij(line_mask) .* (1 - pij(line_mask)));
        
        %Workout which directed pixels are lines
        %theta_line_mask = line_mask(theta_mask);

        %Workout dominant angle in orientation field
        if args.use_lines
            cs = complex(...
                mean(cos(2*ori_map(line_mask))),...
                mean(sin(2*ori_map(line_mask))));
            
            [mu_hat, kappa_hat] = von_mises_mle(2*ori_map(line_mask));
        else
            cs = complex(...
                mean(cos(2*ori_map(r_mask & mam_mask))),...
                mean(sin(2*ori_map(r_mask & mam_mask))));
            [mu_hat, kappa_hat] = von_mises_mle(2*ori_map(r_mask & mam_mask));
        end
        
        %Compute map of Von M probabilities
        %pvm_ij = zeros(800);
        if args.do_vm
            idx = find(line_mask);
            pvm_ij = zeros(length(idx),1);
            for jj = 1:length(idx)
                kk = idx(jj);
                alpha = pvm(kk);
                %pvm_ij(kk) = von_mises_cdf([-alpha alpha], mu_hat - 2*theta(kk), kappa_hat);
                pvm_ij(jj) = von_mises_cdf([-alpha alpha], mu_hat - 2*theta(kk), kappa_hat);
            end

            %Compute p assuming von mises
            %var_out.p_line_vm(ii,aa) = sum(pvm_ij(line_mask)) / var_out.N_line(ii,aa);
            var_out.p_line_vm(ii,aa) = sum(pvm_ij) / var_out.N_line(ii,aa);

            %figure; imagesc(pvm_ij); axis image; colorbar; colorbar;
        end

        %Get angle differences between dominant angle and the angles of
        %pixels directed to centre
        %angle_diff = angle(conj(ori_pm(theta_mask).^2) .* cs) / 2;
        angle_diff = angle(conj(ori_pm.^2) .* cs) / 2;
        
        %Now histogram these differences
        %var_out.n_hist_all(ii, :, aa) = hist(angle_diff, bins);
        var_out.n_hist_all(ii, :, aa) = hist(angle_diff(theta_mask), bins);
        var_out.N_hist_line(ii, :, aa) = hist(angle_diff(line_mask), bins);
        
        %Also histogram the difference for just line points
        %var_out.n_hist_line(ii, :, aa) = hist(angle_diff(theta_line_mask), bins);
        var_out.n_hist_line(ii, :, aa) = hist(angle_diff(line_mask & theta_mask), bins);
        var_out.N_hist_all(ii, :, aa) = hist(angle_diff(r_mask & mam_mask), bins);
        
        %var_out.x_line(ii, :, aa) = sum(1 - pij(line_mask & theta_mask)) - sum(pij(line_mask & theta_mask));
        
        if ~isempty(args.save_path)
            save(args.save_path, 'var_out');
        end
        
    end   
end

var_out.n_all = [sum(var_out.n_hist_all(:,:,1),2) sum(var_out.n_hist_all(:,:,2),2)];
var_out.n_line = [sum(var_out.n_hist_line(:,:,1),2) sum(var_out.n_hist_line(:,:,2),2)];

var_out.f1_all = (var_out.n_all - (var_out.N_all * var_out.p_all)) ./...
    sqrt(var_out.N_all * var_out.p_all * (1 - var_out.p_all));
var_out.f1_line = (var_out.n_line - (var_out.N_line .* var_out.p_line)) ./...
    sqrt(var_out.N_line .* var_out.p_line .* (1 - var_out.p_line));

if args.do_vm
    var_out.f1_line_vm = (var_out.n_line - (var_out.N_line .* var_out.p_line_vm)) ./...
        sqrt(var_out.N_line .* var_out.p_line_vm .* (1 - var_out.p_line_vm));
end