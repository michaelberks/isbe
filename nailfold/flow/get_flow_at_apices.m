function [mean_flow, median_flow, prctile_flow, flow_measurable] = ...
    get_flow_at_apices(vessel_data, mosaic_flow)
%GET_FLOW_AT_APICES *Insert a one line summary here*
%   [mean_flow, median_flow, prctile_flow, flow_measurable] = get_flow_at_apices(vessel_data, mosaic_flow)
%
% Inputs:
%      vessel_data - *Insert description of input variable here*
%
%      mosaic_flow - *Insert description of input variable here*
%
%
% Outputs:
%      mean_flow - *Insert description of input variable here*
%
%      median_flow - *Insert description of input variable here*
%
%      prctile_flow - *Insert description of input variable here*
%
%      flow_measurable - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 06-Apr-2016
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester
if ~isempty(vessel_data.apex_measures.distal)
    distal_xy = vessel_data.apex_measures.distal.apex_xy / ...
        vessel_data.resize_factor;
    distal_xy(:,1) = distal_xy(:,1) + mosaic_flow.padding(1);
    distal_xy(:,2) = distal_xy(:,2) + mosaic_flow.padding(3);
    num_vessels = size(distal_xy,1);
    flow_measurable = false(num_vessels,1);
    mean_flow = zeros(num_vessels,1);
    median_flow = zeros(num_vessels,1);
    prctile_flow = zeros(num_vessels,1);

    for i_ve = 1:num_vessels
        ax = distal_xy(i_ve,:);
        if isnan(ax(1))
            continue;
        end
        apex_mask = bwselect(mosaic_flow.mask, ax(1), ax(2));

        if any(apex_mask(:))
            flow_rates = abs(mosaic_flow.max(apex_mask));
            mean_flow(i_ve) = mean(flow_rates);
            median_flow(i_ve) = median(flow_rates);
            flow_measurable(i_ve) = 1;
            prctile_flow(i_ve) = prctile(flow_rates,90);
        end
    end
end
