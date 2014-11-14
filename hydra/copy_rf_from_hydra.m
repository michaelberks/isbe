function [] = copy_rf_from_hydra(job_id)
%COPY_RF_FROM_HYDRA Copy a random forest constructed and combined on hydra
%from isbe-san1 to local machine hard disk
%   [] = copy_rf_from_hydra(job_id)
%
% Inputs:
%      job_id - *Insert description of input variable here*
%
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 13-Aug-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester 

mb_combine_rfs(...
    'rf_dir', ['Z:\asymmetry_project\data\line_detection_rfs\' num2str(job_id) '\'],...
    'tree_dir',  [num2str(job_id) '\'],...
    'tree_root', [asymmetryroot 'data/line_detection_rfs/'],...
    'replace_tree_root', 'Z:\asymmetry_project\data\line_detection_rfs\',...
    'copy_trees', 1,...% the optional arguments
    'delete_trees', 0,...
    'save_path', [asymmetryroot 'data\line_detection_rfs\' num2str(job_id) '\']);