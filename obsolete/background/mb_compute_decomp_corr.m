function [decomp_corr] = mb_compute_decomp_corr(decomp, decomp_type, varargin)


args = u_packargs(varargin, 0, ...
			'ComputeDecompLocal', 11, ...
			'ComputeDecompSubbands', 1, ...
            'SampleMask', [],...
            'Overlap', 32,...
            'Threshold', 0.01 ...
			);

if isempty(args.SampleMask)
    if strcmp(decomp_type, 'pyr')
        args.SampleMask = ones(size(decomp{2,1}));
    else
        args.SampleMask = ones(size(decomp{1}(:,:,1)));
    end
end

args.SampleMask(1:args.Overlap, :) = 0;
args.SampleMask(end-args.Overlap+1:end, :) = 0;
args.SampleMask(:,1:args.Overlap) = 0;
args.SampleMask(:,end-args.Overlap+1:end) = 0;

        
if args.ComputeDecompLocal
    decomp_corr.local =...
        compute_local_corr(decomp, args.SampleMask, args.Threshold, decomp_type, args.ComputeDecompLocal);
end
if args.ComputeDecompSubbands
    decomp_corr.subbands =...
        compute_subbands_corr(decomp, args.SampleMask, args.Threshold, decomp_type);
end
                

                
function local_corr = compute_local_corr(decomp, sample_mask, thresh, decomp_type, win_size)

if strcmp(decomp_type, 'pyr')
    num_levels = size(decomp,1)-2;
    num_oris = size(decomp,2);
    dim = num_levels*num_oris;
else
    num_levels = size(decomp,1)-1;
    num_oris = size(decomp{1},3); %always 6
    dim = num_levels*num_oris;
end

local_corr.p_val = zeros(dim, win_size^2);
local_corr.rho = zeros(dim, win_size^2);
local_corr.map = zeros(dim, win_size^2);

curr_idx = 0;

[r1_pts c1_pts] = find(sample_mask);

for level = 1:num_levels
    
    r_pts = ceil(r1_pts / 2^(level-1));
    c_pts = ceil(c1_pts / 2^(level-1));
    
    [u_pts] = unique([r_pts c_pts], 'rows');
    
    for ori = 1:num_oris
        
        if strcmp(decomp_type, 'pyr')
            subband = decomp{level+1,ori};
        else
            subband = real(decomp{level}(:,:,ori));
        end

        curr_idx = curr_idx + 1;
        
        data = zeros(size(r_pts,1), win_size^2);
    
        for curr_point = 1:size(u_pts, 1)

            r = u_pts(curr_point, 1);
            c = u_pts(curr_point, 2);
            data(curr_point,:) = ...
                reshape(sample_window(subband, win_size, r, c), 1, []);

        end
            
        %
        if any(isnan(data(:)))
            %Get rid of rows with Nan points (caused by sampling off edge
            %of subband image
            display('overlap too small, removing points')
            [nan_r nan_c] = find(isnan(data)); %#ok
            data(nan_r,:) = [];
        end
        
        [rho p_val] = corr(data);
        
        local_corr.p_val(curr_idx, :) = p_val(0.5*(1+win_size^2),:);
        local_corr.rho(curr_idx, :) = rho(0.5*(1+win_size^2),:);
        local_corr.map(curr_idx, :) = local_corr.p_val(curr_idx, :) < thresh;
        clear data;
        
    end
end


function subbands_corr = compute_subbands_corr(decomp, sample_mask, thresh, decomp_type)

if strcmp(decomp_type, 'pyr')
    num_levels = size(decomp,1)-2;
    num_oris = size(decomp,2);
    dim = num_levels*num_oris;
else
    num_levels = size(decomp,1)-1;
    num_oris = size(decomp{1},3); %always 6
    dim = 2*num_levels*num_oris;
end

[r1_pts c1_pts] = find(sample_mask);

data = zeros(size(r1_pts,1), dim);

curr_idx = 0;

for level = 1:num_levels 
    
    
    r_pts = ceil(r1_pts/2^(level-1));
    c_pts = ceil(c1_pts/2^(level-1));

    if strcmp(decomp_type, 'pyr')
        idx = sub2ind(size(decomp{level+1,1}), r_pts, c_pts);
    else
        idx = sub2ind(size(decomp{level}(:,:,1)), r_pts, c_pts);
    end
    
    for ori = 1:num_oris
        
        curr_idx = curr_idx + 1;
        if strcmp(decomp_type, 'pyr')
            subband = decomp{level+1,ori};
            data(:,curr_idx) = subband(idx);         
        else
            subband = decomp{level}(:,:,ori);
            data(:,curr_idx) = real(subband(idx));
            data(:,curr_idx+num_levels*num_oris) = imag(subband(idx));
        end
   
    end 
end

clear r* c* sample_mask idx*


%
%[subbands_corr.rho subbands_corr.p_val] = corr(data);
cov_mat = cov(data);
sd_mat = sqrt(diag(cov_mat)) * ones(1, size(data,2));
subbands_corr.rho = cov_mat ./ (sd_mat .* sd_mat');
%subbands_corr.map = subbands_corr.p_val < thresh;
clear data;