function [smooth_bg] = mb_smooth_mammo_bg(mammo_bg,levels)
%MB_SMOOTH_MAMMO_BG *Insert a one line summary here*
%   [smooth_bg] = mb_smooth_mammo_bg(mammo_bg,levels)
%
% Inputs:
%      mammo_bg- *Insert description of input variable here*
%
%      levels- *Insert description of input variable here*
%
%
% Outputs:
%      smooth_bg- *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 24-Feb-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

mammo_dt = dtwavexfm2(mammo_bg, levels+1);
%%
mammo_dt2 = mammo_dt;
mammo_dt3 = mammo_dt;

for lev = 1:levels
    for ori = 1:6
        subband = mammo_dt2{lev}(:,:,ori);
        [dummy sort_idx] = sort(abs(subband(:)));
        sort_idx(end-round(2*end/5):end) = [];
        subband(sort_idx) = 0;
        mammo_dt2{lev}(:,:,ori) = subband;
    end
    mammo_dt3{lev}(:) = 0;
end
mammo_bg2 = dtwaveifm2(mammo_dt2);

mammo_noise = mammo_bg - mammo_bg2;
mammo_coarse = dtwaveifm2(mammo_dt3);

smooth_bg = mammo_coarse + mammo_noise;
