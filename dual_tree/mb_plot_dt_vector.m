function [] = mb_plot_dt_vector(dt_vector)
%MB_PLOT_DT_VECTOR *Insert a one line summary here*
%   [] = mb_plot_dt_vector(dt_vector)
%
% Inputs:
%      dt_vector- *Insert description of input variable here*
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
% Created: 03-Nov-2008
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester 

%This daft bit of code seems to be only way of generating a suitably
%large circular grid!
h = polar([0 2*pi], [0 max(dt_vector(1:3:end-2))]); delete(h);


levels = length(dt_vector) / 3;

colors = 'rgbmcky';

hold on;
title('Polar plot of the amplitude and phase at each scale');

legend_labels = cell(1,2*levels);
for lev = 1:levels
    
    compass(dt_vector(lev*3-2)*exp(i*(dt_vector(lev*3-1)+dt_vector(end))), colors(lev));
    compass(dt_vector(lev*3-2)*exp(i*dt_vector(lev*3-1)), [colors(lev), ':']);
    legend_labels{2*lev - 1} = ['Level ', num2str(lev)];
    legend_labels{2*lev} = ['Level ', num2str(lev), ' relative to theta'];
    
end
legend(legend_labels);