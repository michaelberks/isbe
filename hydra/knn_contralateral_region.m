function knn_contralateral_region(part_idx, num_parts)

%assume num_parts is square
num_parts = sqrt(num_parts);

contra_pair = u_load([asymmetryroot 'data/contralateral/024RML.mat']);
map_args.region1 = imresize(contra_pair.abnormal_roi, 0.5, 'bilinear');
map_args.region2 = imresize(contra_pair.normal_roi, 0.5, 'bilinear');
clear contra_pair;

[row col] = size(map_args.region1);
[row_idx col_idx] = ind2sub([num_parts num_parts], part_idx); 

sr = round((row_idx - 1) * row / num_parts) + 1;
sc = round((col_idx - 1) * col / num_parts) + 1;

er = round(row_idx * row / num_parts);
ec = round(col_idx * col / num_parts);

map_args.mask = false(row, col);
map_args.mask(sr:er, sc:ec) = 1;

map_args.num_samples2 = 1e4;
map_args.k = 10;
map_args.win_size = 3;
map_args.num_levels = 4;

dist_map = knn_2_region_map(map_args);
save([asymmetryroot 'results/knn_maps/024RML_knn_map_W3L5N1e4K10_' zerostr(part_idx, floor(2*log10(num_parts))+1) '.mat'], 'dist_map');