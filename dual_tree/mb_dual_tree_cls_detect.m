function [cls] = mb_dual_tree_cls_detect(mam_region, num_levels, cutoffs, percents, phase_thresh, if_plot)

if nargin < 6
    if_plot = 0;
end
if nargin < 5
    phase_thresh = [-inf -inf];
end

dual_tree = dtwavexfm2(mam_region, num_levels+1);
[dt_ilp dt_icp] = mb_dual_tree_transform(dual_tree);

cls.map = cell(num_levels,1);
cls.orientations = cell(num_levels,1);

for level = 1:num_levels
    
   cls.orientations{level} = compute_max_icp(dt_icp{level});
   icp_band = max(dt_ilp{level}, [], 3);
   cls.map{level} =...
       icp_hysterisis(cls.orientations{level}, icp_band, [cutoffs(level,percents(1)) cutoffs(level,percents(2))], phase_thresh, if_plot);
    
end