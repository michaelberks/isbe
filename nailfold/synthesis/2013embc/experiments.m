clear;

n_frames = 180;

sigma_x = []; sigma_y = []; sigma_w = [];

% clean sequence
bg_weight = 0;

% First element always default
flowrate_seq = [1.0, 0.5:0.25:2.0]; flowrate = 1.0;
n_cells_seq = round([640, 40*2.^(0:0.25:5)]); n_cells = 640;
contrast_seq = [8:8:32,48,64]; contrast = 24;
j_sigma_seq = 0.0:0.5:2.5; j_sigma = 0.0;
b_sigma_seq = [0.0:0.2:1.0, 1.5:0.5:4.0]; b_sigma = 0.0;
c_sigma_seq = 0.0:0.025:0.25; c_sigma = 0.0;
sigma_n_seq = 0.0:0.5:2.0; sigma_n = 0.0;

flowrate_seq = [0.5, 1.0, 2.0]; flowrate = 1.0;
n_cells_seq = round([320, 640, 1280]); n_cells = 640;
b_sigma_seq = [0, 1, 2]; b_sigma = 0.0;
c_sigma_seq = [0.0, 0.1, 0.25]; c_sigma = 0.0;

% noisy sequence
% j_sigma = 2.0; b_sigma = 3.0; c_sigma = 0.1; sigma_n = 2.0; bg_weight = 20;

% flowsubdir = 'flowrate'; for flowrate = flowrate_seq
% flowsubdir = 'n_cells'; for n_cells = n_cells_seq
% flowsubdir = 'contrast0'; for contrast = 0 % contrast_seq
% flowsubdir = 'jitter_sigma'; for j_sigma = j_sigma_seq
% flowsubdir = 'brightness_sigma'; for b_sigma = b_sigma_seq(2:end)
flowsubdir = 'contrast_sigma'; for c_sigma = c_sigma_seq(2:end)
% flowsubdir = 'sigma_n'; for sigma_n = sigma_n_seq([end-2])
    synthesize_images('n_frames', n_frames, ...
                      'n_cells', n_cells, ...
                      'sigma_x', sigma_x, ...
                      'sigma_y', sigma_y, ...
                      'sigma_w', sigma_w, ...
                      'max_flow', flowrate, ...
                      'jitter_sigma', j_sigma, ...
                      'brightness_sigma', b_sigma, ...
                      'contrast_sigma', c_sigma, ...
                      'bg_weight', bg_weight, ...
                      'sigma_n', sigma_n, ...
                      ...
                      'brightness0', 160, ...
                      'contrast0', contrast);
%     pt_flow_sandpit;
end
