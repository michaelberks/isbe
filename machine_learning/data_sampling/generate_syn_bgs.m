function [] = generate_syn_bgs(bg_dir, syn_dir, varargin)
%GENERATE_SYN_BGS *Insert a one line summary here*
%   [] = generate_syn_bgs(varargin)
%
% GENERATE_SYN_BGS uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%
% Optional Arguments:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 27-Jan-2011
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin, '0', ...
    'image_list', [],...
    'image_fmt', '.mat',...
    'num_levels', 6,...
    'discard_pcnt', 0.4,...
    'plot', 0);
clear varargin;

%check dirs are filesep terminated
if ~strcmp(bg_dir(end), '/')
    bg_dir(end+1) = '/';
end
%check dirs are filesep terminated
if ~strcmp(syn_dir(end), '/')
    syn_dir(end+1) = '/';
end

%get list of images if not supplied by user
if isempty(args.image_list)
    args.image_list = dir([bg_dir '*' args.image_fmt]);
end

%make output dir if it doesn't exist
if ~exist(syn_dir, 'dir')
    mkdir(syn_dir);
end

%Do strcmp now to usual logical compare in loop
do_mat = strcmp(args.image_fmt, '.mat');
%Loop through the real backgrounds
for ii = 1:length(args.image_list)

    %load the real background
    if do_mat
        real_bg = double(u_load([bg_dir args.image_list(ii).name]));
    else
        real_bg = double(imread([bg_dir args.image_list(ii).name]));
    end
    
    %Compute its dual tree decomp
    real_dt = dtwavexfm2(real_bg, args.num_levels+1);
    
    %Make copy of the dual-tree
    real_dt2 = real_dt;

    %Loop throught the required levels
    for lev = 1:args.num_levels
        for ori = 1:6
            
            %Extract each oriented sub-bands from the first copy
            subband = real_dt{lev}(:,:,ori);
            
            %Sort the wavelet coefficients by magnitude
            [dummy sort_idx] = sort(abs(subband(:)));
            
            %Discard the lowest X percent by setting to zero
            sort_idx(end-round(args.discard_pcnt*end):end) = [];
            subband(sort_idx) = 0;
            
            %Replace the modified sub-band
            real_dt{lev}(:,:,ori) = subband;
        end
        
        %In the second copy set all coefficients to zero
        real_dt2{lev}(:) = 0;
    end

    %Reconstruct the bg from the 1st copy - this will have noise smoothed
    %out
    real_smooth = dtwaveifm2(real_dt);

    %Subtract the smooth bg from the original to get a map of noise
    real_noise = real_bg - real_smooth;
    
    %Reconstruct from the 2nd dual-tree copy to get the low-level variation
    %in the bg
    real_coarse = dtwaveifm2(real_dt2);

    %Add the noise map to the coarse bg to obtain the synthetic background
    syn_bg = real_coarse + real_noise;
    
    if args.plot && ii < 20
        figure; 
        subplot(2,2,1); imagesc(real_smooth); axis image; colormap(gray(256));
        subplot(2,2,2); imagesc(real_noise); axis image; colormap(gray(256));
        subplot(2,2,3); imagesc(real_coarse); axis image; colormap(gray(256));
        subplot(2,2,4); imagesc(syn_bg); axis image; colormap(gray(256));
    end
     
    %save the synthetic bg
    if do_mat
        save([syn_dir 'bg', zerostr(ii,5) args.image_fmt], 'syn_bg');
    else
        imwrite(syn_bg, [syn_dir 'bg', zerostr(ii,5) args.image_fmt]);
    end
    
end