nailfold_dir = 'I:\nailfold\software\ncm markup tool\ncm_markup\Inbox\';
nailfold_list = dir([nailfold_dir '*mosaic.bmp']);
num_n = length(nailfold_list);
%%
mkdir([nailfold_dir 'g2d\analytic\orientations\']);
mkdir([nailfold_dir 'g2d\analytic\responses\']);
mkdir([nailfold_dir 'g2d\analytic\scales\']);

for ii = 1:length(nailfold_list)
    
    nailfold = imresize(imread([nailfold_dir nailfold_list(ii).name]), 0.5, 'bilinear');
    [response_map ori_map scale_map] = gaussian_2nd_derivative_line(nailfold, [1 2]);
    save_uint8([nailfold_dir 'g2d\analytic\orientations\' nailfold_list(ii).name(1:end-4) '_ori.mat'], ori_map);
    save([nailfold_dir 'g2d\analytic\responses\' nailfold_list(ii).name(1:end-4) '_response.mat'], 'response_map');
    save_uint8([nailfold_dir 'g2d\analytic\scales\' nailfold_list(ii).name(1:end-4) '_scale.mat'], scale_map);
end
%%
mkdir([nailfold_dir 'masks']);
for ii = 1:num_n
    nailfold = imresize(imread([nailfold_dir nailfold_list(ii).name]), 0.5, 'bilinear');
    edge_mask = nailfold >= 254;
    edge_label = bwlabel(edge_mask, 4);
    
    counts = hist(edge_label(:), 0:max(edge_label(:)));
    [~, lab_idx] = max(counts(2:end));
    edge_mask = edge_label == lab_idx;
    
    [vessels apex_data] = read_vessels_to_excel([nailfold_dir nailfold_list(ii).name(1:end-4) '_markup.txt']);
    apex_data = apex_data/2;
    remove_rows = ~apex_data(:,5) | ~apex_data(:,6);
    apex_data(remove_rows,:) = [];
    
    if isempty(apex_data); 
        region_mask = [];    
    else 
        p1 = polyfit(apex_data(:,5), apex_data(:,6), 1);
        p2 = polyfit(apex_data(:,5), apex_data(:,6), 2);
        x = 1:size(nailfold, 2);
        y1 = polyval(p1, x);
        y2 = polyval(p2, x);
        
        y_diff = polyval(p2, apex_data(:,5)) - apex_data(:,6);
        y_off = max(max(y_diff), 0);
        y2u = y2 - y_off;
        y2l = y2 + 75;

        region_mask = poly2mask([x x(end:-1:1)], [y2u y2l(end:-1:1)], size(nailfold,1), size(nailfold,2));
        region_mask = region_mask & ~edge_mask;

        nailfold_masked = nailfold;
        %nailfold_masked(~region_mask) = mean(nailfold(region_mask));

        if ii < 21
            figure; 
            subplot(2,1,1); imgray(nailfold_masked); 
            %plot(x, y1, 'g');
            plot(x, y2, 'y--');
            plot(x, y2u, 'y-');
            plot(x, y2l, 'y-');
            plot(apex_data(:,5),apex_data(:,6), 'rx');
            subplot(2,1,2); imgray(region_mask); 
            plot([x x(end:-1:1)], [y2u y2l(end:-1:1)]);
        end
    end
    
    save([nailfold_dir 'masks/' nailfold_list(ii).name(1:end-4) '_mask.mat'], 'region_mask');
end
%%
r_vals = zeros(num_n*1000,1);
for ii = 1:num_n
    nailfold = imresize(imread([nailfold_dir nailfold_list(ii).name]), 0.5, 'bilinear');
    response_map = u_load([nailfold_dir 'g2d\analytic\responses\' nailfold_list(ii).name(1:end-4) '_response.mat']);
    mask = u_load([nailfold_dir 'masks/' nailfold_list(ii).name(1:end-4) '_mask.mat']);
    %response_map(response_map < 0) = 0;
    
    if isempty(mask)
        r_vals(offset+(1:1000)) = NaN;
    else
        mask = mask & (response_map > 0);
        responses = response_map(mask);
        shuffle = randperm(numel(responses));
        offset = (ii-1)*1000;
        r_vals(offset+1) = min(responses(:));
        r_vals(offset+2) = max(responses(:));
        r_vals(offset+(3:1000)) = responses(shuffle(1:998));
    end
end
r_vals(isnan(r_vals)) = [];
%
figure; hist(r_vals, 1000);
r_vals = unique(r_vals);
probs = linspace(0,1,numel(r_vals));
figure; plot(r_vals, probs);
%%
mkdir([nailfold_dir 'vessel_probs']);
for ii = 1:num_n
    nailfold = imresize(imread([nailfold_dir nailfold_list(ii).name]), 0.5, 'bilinear');
    response_map = u_load([nailfold_dir 'g2d\analytic\responses\' nailfold_list(ii).name(1:end-4) '_response.mat']);
    mask = u_load([nailfold_dir 'masks/' nailfold_list(ii).name(1:end-4) '_mask.mat']);
    
    if isempty(mask); continue; end
    
    response_map_p = interp1(r_vals, probs, response_map, 'cubic');
    mask = u_load([nailfold_dir 'masks/' nailfold_list(ii).name(1:end-4) '_mask.mat']);
    response_map_p(response_map_p < 0) = 0;
    response_map_p(response_map_p > 1) = 1;
    response_map_p(response_map < 0) = 0;
    response_map_p(~mask) = 0;
    
    if ii < 21
        figure;
        subplot(2,1,1); imgray(nailfold);
        subplot(2,1,2); imgray(response_map_p);
    end
    save([nailfold_dir 'vessel_probs\' nailfold_list(ii).name(1:end-4) '_vp.mat'], 'response_map_p');
end
%%
entropies = zeros(num_n,1);
mkdir([nailfold_dir 'entropy']);
for ii = 1:num_n
    [vessels apex_data] = read_vessels_to_excel([nailfold_dir nailfold_list(ii).name(1:end-4) '_markup.txt']);
    apex_data = apex_data/2;
    remove_rows = ~apex_data(:,5) | ~apex_data(:,6);
    apex_data(remove_rows,:) = [];
    
    if isempty(apex_data); continue; end
    
    nailfold = imresize(imread([nailfold_dir nailfold_list(ii).name]), 0.5, 'bilinear');
    response_map_r = u_load([nailfold_dir 'vessel_probs\' nailfold_list(ii).name(1:end-4) '_vp.mat']);
    response_map = u_load([nailfold_dir 'g2d\analytic\responses\' nailfold_list(ii).name(1:end-4) '_response.mat']);
    orientation_map = load_uint8([nailfold_dir 'g2d\analytic\orientations\' nailfold_list(ii).name(1:end-4) '_ori.mat']);
    
    vessel_mask = bwselect(response_map_r > 0.1, apex_data(:,5), apex_data(:,6), 8);
%     vessel_nms = bwareaopen(...
%             mb_non_maximal_supp(response_map, orientation_map) > 0, 10);
%     [y_vessels x_vessels] = find(vessel_nms);    
%     ori_all = mod(floor(180*orientation_map(vessel_mask)/pi),180) + 1;
%     res_all = response_map_r(vessel_mask);
%     ori_hist = full( sparse(ori_roi,1,res_roi,180,1) );
    ori_hist = hist(mod(orientation_map(vessel_mask),pi), linspace(0,pi,90));
    total_entropy = mb_entropy(ori_hist);
    entropies(ii) = total_entropy;
    
    if ii < 21
        
        response_map_r(~vessel_mask) = 0;
        combined_map = response_map_r .* exp(2*1i*orientation_map);
        figure; 
        a1 = subplot(2,4,1:3); imgray(nailfold); plot(apex_data(:,5),apex_data(:,6), 'rx');
        a2 = subplot(2,4,5:7); imgray(complex2rgb(combined_map));
        subplot(2,4,[4 8]); plot(linspace(0,pi,90), ori_hist); title(['Entropy: ' num2str(total_entropy)] );
        imwrite(complex2rgb(combined_map), [nailfold_dir 'entropy/' nailfold_list(ii).name(1:end-4) '_weighted_ori.bmp']);
    end
        
        
%     save([nailfold_dir 'entropy/' nailfold_list(ii).name(1:end-4) '_entropy.mat'], 'entropy');
%     
%     entropies = zeros(size(apex_data,1),1);
%     
%     for jj = 1:size(apex_data,1)
%         x = apex_data(jj,5);
%         y = apex_data(jj,6);
%         x1 = x - 40;
%         x2 = x + 40;
%         y1 = y - 25;
%         y2 = y + 75;
%         
%         
% %         [~, idx] = min((x_vessels-x).^2 + (y_vessels-y).^2);
% %         this_vessel = bwselect(vessel_nms, x_vessels(idx), y_vessels(idx), 8);
% %         [yv xv] = find(this_vessel);
%         
%         plot(a1, [x1 x2 x2 x1 x1], [y1 y1 y2 y2 y1], 'g');
%         plot(a2, [x1 x2 x2 x1 x1], [y1 y1 y2 y2 y1], 'w', 'linewidth', 2);
% %         plot(a1, xv, yv, 'r.');
%         
%         roi = poly2mask([x1 x2 x2 x1], [y1 y1 y2 y2], size(nailfold,1), size(nailfold,2));
%         
%         ori_roi = mod(floor(180*orientation_map(roi)/pi),180) + 1;
%         res_roi = response_map_r(roi);
%         %ori_hist = full( sparse(ori_roi,1,res_roi,180,1) );
%         ori_hist = hist(mod(orientation_map(roi),pi), linspace(0,pi,90));
%         entropies(jj) = mb_entropy(ori_hist);
%         %figure; bar((0:179)', ori_hist);       
%     end
%     title(a1, ['Mean entropy: ' num2str(mean(entropies))] );
    
end
%%
names = cell(num_n,1);
for ii = 1:num_n
    names{ii} = nailfold_list(ii).name(1:end-4);
end
xlswrite('test.xlsx', names, 1, 'A2');
xlswrite('test.xlsx', entropies, 1, 'B2');
%%
show_ori_wheel

    