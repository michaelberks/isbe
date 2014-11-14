function quick_apex_marking(contour_dir, vessel_selection)



vessel_list = dir([contour_dir, '*contour.mat']);


if nargin < 2
    num_vessels = length(vessel_list);
    vessel_selection = 1:num_vessels;
else
    num_vessels = length(vessel_selection);
end

screen_size = get(0,'ScreenSize');
seg_fig = figure(...
    'windowstyle', 'normal',...
    'position', screen_size);

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
    
    contour_struc = load([contour_dir vessel_list(i_ve).name]);
    
    if isfield(contour_struc, 'apex_idx')
        selection = questdlg('Markup already exists. Process image anyway?',...
             'Warning',...
             'Yes','No','Yes');           
        if strcmpi(selection, 'No')
            continue;
        end
    end    
    
    outer_edge = contour_struc.outer_edge;
    inner_edge = contour_struc.inner_edge;
    vessel_centre = contour_struc.vessel_centre;
    
    %Display contour
    figure(seg_fig);
    axis equal ij; hold all;
    set(seg_fig, 'Name', ['Vessel ', num2str(counter), ' of ', num2str(num_vessels),...
        ': ', vessel_list(i_ve).name]);
    
    plot(outer_edge(:,1), outer_edge(:,2), 'g');
    plot(inner_edge(:,1), inner_edge(:,2), 'b');
    plot(vessel_centre(:,1), vessel_centre(:,2), 'k');
    plot(...
            [outer_edge(:,1) inner_edge(:,1)]',...
            [outer_edge(:,2) inner_edge(:,2)]', 'y');
        
    xlabel('Click on apex');
    
    happy = false;
    while ~happy
        %Get user to click on apex
        [x, y] = ginput();
        num_pts = length(x);
        
        apex_idx = zeros(num_pts,1);
        
        for i_pt = 1:num_pts
            %Find nearest vessel centre point to click
            dists = (vessel_centre(:,1) - x(i_pt)).^2 + (vessel_centre(:,2) - y(i_pt)).^2;
            [~, apex_idx(i_pt)] = min(dists);
        end
        plot(vessel_centre(apex_idx,1), vessel_centre(apex_idx,2), 'rx');
        
        answer = questdlg('Are the points ok?',...
             'Happy?',...
             'Yes','No','Yes');  
        happy = strcmpi(answer, 'Yes');
    end
    
    
    save([contour_dir vessel_list(i_ve).name],...
        'vessel_centre', 'outer_edge', 'inner_edge', 'apex_idx');          
end

close(seg_fig);
display(['Manual segmentation terminated at image ', num2str(i_ve)]);
            