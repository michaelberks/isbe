function [cls_map] =  mb_cls_map_orientations(cls_result, orientations)

% Create CLS map structure from CLS_results structure (as returned from
% MB_CLS_SELECTION), so that each level of the CLS map structure defines
% pixels marked as CLS in that OR any coarser level

if nargin < 2
    orientations = [2*pi/5 -2*pi/5; pi/5 2*pi/5; 0 pi/5; -pi/5 0; -2*pi/5 -pi/5];
end

levels = size(cls_result, 1);
cls_map = cell(levels, 1);
num_orientations = size(orientations, 1);

for curr_level = 1:levels
    
    %pre-allocate map array for this level
    cls_map{curr_level} = zeros(size(cls_result{curr_level}.CLS));
    
    %Quantise the orientations
    
    %First orientation is assumed to overlap cyclical range of rads => '|'
    ori_map = cls_result{curr_level}.CLS & ...
                (cls_result{curr_level}.GaborOrientation >= orientations(1, 1) | ...
                cls_result{curr_level}.GaborOrientation < orientations(1, 2));
    cls_map{curr_level}(ori_map) = 1;
    
    %All other orientations need '&'
    for ori = 2:num_orientations
        ori_map = cls_result{curr_level}.CLS & ...
                    cls_result{curr_level}.GaborOrientation >= orientations(ori, 1) & ...
                    cls_result{curr_level}.GaborOrientation < orientations(ori, 2);
        
        cls_map{curr_level}(ori_map) = ori;       
    end
        
end

