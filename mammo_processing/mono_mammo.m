function [] = mono_mammo(im_in, no_filt, min_wav, mult, sigma)
% Applies the monogenic signal to a mammograms and plots images of the
% phase and orientation maps at multiple scales
%
%notes: quick test function exploring the monogenic signal processing described
%Pan et al. as a method for segmenting mammograms

[local_amp, local_ori, local_phase] = monogenic(im_in, no_filt, min_wav, mult, sigma, 1);
[rows cols] = size(im_in);

for ii = 1:4
    
    phase_mask = local_phase{ii} > pi/2;
    ori_map = zeros(rows, cols);
    ori_map(local_ori{ii} > -pi/2) = 1;
    ori_map(local_ori{ii} > 0) = 2;
    ori_map(local_ori{ii} > pi/2) = 3;
    
    subplot(4, 4, ii); imagesc(local_phase{ii}); axis image; colorbar;
    subplot(4, 4, 4+ii); imagesc(phase_mask); axis image; colorbar;
    subplot(4, 4, 8+ii); imagesc(ori_map); axis image; colorbar;
    
    clear ori_map;
    clear phase_mask;
end