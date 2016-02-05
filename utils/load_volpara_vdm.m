function [norm_vdm, H, raw_vdm] = load_volpara_vdm(filepath, normalise_by_width)
%LOAD_VOLPARA_VDM *Insert a one line summary here*
%   [norm_vdm, H, raw_vdm] = load_volpara_vdm(filepath, normalise_by_width)
%
% Inputs:
%      filepath - path to Volpara VDM
%
%      normalise_by_width - {0,1} if true, get height from filename and
%      normalise using height
%
%
% Outputs:
%      norm_vdm - normalised 
%
%      H - breast thickness, as recorded in file name
%
%      raw_vdm - original unormalised VDM
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 09-Nov-2015
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester

if exist('normalise_by_width', 'var') && normalise_by_width
    [~, fname, ~] = fileparts(filepath);
    s_idx = strfind(fname, 'Map_H_') + 6;
    e_idx = strfind(fname, 'mm') - 1;
    
    if isempty(s_idx) || isempty(e_idx)
        warning(['Breast thickness not found in ' fname '. Cannot normalise by thickness']);
        H = 1;
    else
        H = str2double(fname(s_idx:e_idx));
    end
else
    H = 1;
end

raw_vdm = imread(filepath);
norm_vdm = (double(raw_vdm)./ 65535 - 0.5) .* 2.0 * H;
norm_vdm(norm_vdm < 0) = 0;
