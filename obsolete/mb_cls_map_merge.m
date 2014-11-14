function [cls_map] =  mb_cls_map_merge(cls_result)

% Create CLS map structure from CLS_results structure (as returned from
% MB_CLS_SELECTION), so that each level of the CLS map structure defines
% pixels marked as CLS in that OR any coarser level

levels = size(cls_result, 1);
cls_map = cell(levels, 1);

for curr_level = levels:-1:1
    
    cls_map{curr_level} = cls_result{curr_level}.CLS;
    
    [r c] = size(cls_map{curr_level});
    
    for lower_level = levels:-1:curr_level+1
        cls_map{curr_level} = cls_map{curr_level} | ...
            imresize(cls_result{lower_level}.CLS, [r c]);
    end
end

