function [] = pt_testbench()
%PT_TESTBENCH Testbench for hydra.
%   [] = pt_testbench()
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
% Created: 23-Feb-2011
% Author: Phil Tresadern 
% Email : philip.tresadern@manchester.ac.uk 
% Phone : +44 (0)161 275 5114 
% Copyright: (C) University of Manchester 

clc; clear all; close all hidden;

% use NAG by default
use_nag = 1;

if ispc && strcmp(get_username,'ptresadern')
	use_nag = 0;
	profile on; profile clear;
end

% build a single tree
predictor = build_rf_line_detector(...
	'rand_seed', [], ...
	'sampling_method','generate_training_data',...
	'split_criterion','dabs',...
	'var_criterion','mabs',...
	'image_type','line',...
	'bg_fmt', 'mat',...
	'bg_stem', [],...
	'num_samples', 500, ...
	'pts_per_image', 50, ...
	'detection_type', 'linear_regression', ...
	'decomp_type', 'dt',...
	'feature_type', 'conj', ...
	'win_size', 1,...
	'n_trees', 1,...
	'd', [],...
	'minimise_size', 0, ...
	'split_min', 10,...
	'w_prior', 0.05, ...
	'use_nag', use_nag);

% check that we have a valid job
forest_job = unixenv('JOB_ID',strtok(predictor.tree_dir,'/'));
if isempty(forest_job)
	error('Forest job not found');
end

if	~isfield(predictor,'regression_method') || ...
	strcmp(predictor.regression_method,'random_forest')
	% build a forest from single tree
	combine_hydra_rfs(forest_job,'line_orientation_rfs');
	forest_dir = 'line_orientation_rfs';
else
	forest_dir = 'line_linear_regression_rfs';
end

% test the forest on a limited dataset
classify_image_set(forest_job,'synthetic_lines/philtres/onesample',...
	'forest_dir',forest_dir);

if ispc
	profstat = profile('status');
	if strcmp(profstat.ProfilerStatus,'on')
		profile report; profile off;
	end
end
	


