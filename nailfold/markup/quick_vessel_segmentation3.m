function quick_vessel_segmentation3(vessel_dir, contour_dir, vessel_selection)



vessel_list = dir([vessel_dir, '*vessel.mat']);


if nargin < 3
    num_vessels = length(vessel_list);
    vessel_selection = 1:num_vessels;
else
    num_vessels = length(vessel_selection);
end

screen_size = get(0,'ScreenSize');
seg_fig = figure(...
    'windowstyle', 'normal',...
    'position', screen_size);
buff = 30;

counter = 0;
for i_ve = vessel_selection
    counter = counter+1;
    selection = questdlg('Continue to next image?',...
                     'Warning',...
                     'Yes','No','Yes');           
    if strcmpi(selection, 'No')
        break
    end
    delete(get(seg_fig, 'children'));
    
    vessel_name = [vessel_dir, vessel_list(i_ve).name];
    save_name = [contour_dir vessel_list(i_ve).name(1:end-4) '_edge.mat'];
    
    if exist(save_name, 'file')
        selection = questdlg('Markup already exists. Process image anyway?',...
             'Warning',...
             'Yes','No','Yes');           
        if strcmpi(selection, 'No')
            continue;
        end
    end    
        
    %Load vessel and make contrast equalised patch
    vessel_struc = u_load(vessel_name);
    v_pts = vessel_struc.v_pts;
    vessel_patch = double(vessel_struc.vessel_patch);
    g = gaussian_filters_1d(16, 48);
    g = g / sum(g);
    im_edges = conv2(g', g, ones(size(vessel_patch)), 'same');
    vessel_patch_smoothed = conv2(g', g, vessel_patch, 'same') ./ im_edges;
    vessel_patch_equalised = vessel_patch - vessel_patch_smoothed;

    %Display vessel patch and annotated pts
    figure(seg_fig);
    set(seg_fig, 'Name', ['Vessel ', num2str(counter), ' of ', num2str(num_vessels),...
        ': ', vessel_list(i_ve).name]);
    
    imgray(vessel_patch_equalised(buff+1:end-buff, buff+1:end-buff)');
    caxis(prctile(vessel_patch_equalised(:), [5 95]));
    plot(v_pts(:,2)-buff, v_pts(:,1)-buff, 'y--');
    plot(v_pts(1,2)-buff, v_pts(1,1)-buff, 'rx');
    plot(v_pts(end,2)-buff, v_pts(end,1)-buff, 'gx');
    
    selection = questdlg('Process image?',...
                     'Warning',...
                     'Yes','No','Yes');           
    if strcmpi(selection, 'No')
        vessel_struc.quality = 2;
        save(vessel_name, 'vessel_struc');
        continue;
    end
    
    if vessel_struc.quality == 2
        vessel_struc.quality = 3;
        save(vessel_name, 'vessel_struc');
    end
    
    %Get outer edge
    loop = true;
    while loop
        [~, xi, yi] = roipoly;
        if ~isempty(xi)
            loop = false;
        end
    end
    outer_edge = [xi yi];
    [~, i] = unique(outer_edge, 'rows', 'first');
    outer_edge = outer_edge(sort(i),:);
    [outer_edge] = spline_contour(outer_edge, [], 2);
    
    plot(outer_edge(:,1), outer_edge(:,2), 'g-');
    
    %Get inner edge
    loop = true;
    while loop
        [~, xi, yi] = roipoly;
        if ~isempty(xi)
            loop = false;
        end
    end            
    inner_edge = [xi yi];
    [~, i] = unique(inner_edge, 'rows', 'first');
    inner_edge = inner_edge(sort(i),:);
    [inner_edge] = spline_contour(inner_edge, [], 2);
    
    plot(inner_edge(:,1), inner_edge(:,2), 'g-');
    
    %Convert coords to aoriginal patch frame and save
    outer_edge = [outer_edge(:,2) outer_edge(:,1)]+buff;
    inner_edge = [inner_edge(:,2) inner_edge(:,1)]+buff;
    
    save(save_name, 'outer_edge', 'inner_edge');          
end

close(seg_fig);
display(['Manual segmentation terminated at image ', num2str(i_ve)]);
            