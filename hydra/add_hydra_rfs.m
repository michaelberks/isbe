function [] = add_hydra_rfs(forest_job, forest_number, detection_type)
%ADD_HYDRA_RFS *Insert a one line summary here*
%   [] = add_hydra_rfs()
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

if nargin < 3
    %Choose between line classification and orientation regression
    if ispc        
        detection_type = 'detection';
    else
        [z detection_type] = unix('echo $DETECTION_TYPE'); detection_type(end) = [];
    end
end

forest_number = zerostr(forest_number, 2);
source_rf_path = [asymmetryroot 'data/line_' detection_type '_rfs/' forest_job '/random_forest' forest_number '.mat'];
target_rf_path = [asymmetryroot 'data/line_' detection_type '_rfs/' forest_job '/random_forest.mat'];


add_rfs(source_rf_path, target_rf_path, 1);

display(['RF ' forest_number ' successfully added to ' asymmetryroot 'data/line_' detection_type '_rfs/' forest_job '/']);
