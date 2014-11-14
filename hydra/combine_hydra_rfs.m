function [] = combine_hydra_rfs(varargin)
%COMBINE_HYDRA_RFS *Insert a one line summary here*
%   [] = combine_hydra_rfs()
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
% Created: 12-Aug-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester
args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    'model_id',             unixenv('MODEL_PATH', 'model'),...
    'predictor_name',       unixenv('PREDICTOR_NAME','predictor'), ...
    'model_root',           [unixenv('DATA_ROOT',[]) unixenv('MODEL_ROOT', [asymmetryroot,'data/models/vessel'])],...
    'make_sampled_maps',     unixenv('MAKE_SAMPLED_MAPS', 0));

model_root = args.model_root;
predictor_name = args.predictor_name;

if args.make_sampled_maps
    sampled_maps_dir = [args.model_id '/sampled_maps/'];
else
    sampled_maps_dir = [];
end
combine_rfs(...
    'rf_dir', [model_root '/' args.model_id '/'],...
    'rf_name', predictor_name,...
    'tree_dir', [args.model_id '/trees/'],...
    'tree_root', [model_root '/'],...
    'copy_trees', 1,...% the optional arguments
    'delete_trees', 1,...
    'sampled_maps_dir', sampled_maps_dir, ...
    'save_path', [model_root '/' args.model_id '/' predictor_name '.mat']);

display(['RF ' args.model_id ' successfully combined and saved to ' model_root '/' args.model_id '/' predictor_name '.mat']);
