function quick_fov_markup(vessel_dir, fov_dir, vessel_selection)



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
    save_name = [fov_dir vessel_list(i_ve).name(1:end-4) '_fov_mask.mat'];
    
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
    
    i1 = imgray(vessel_patch_equalised');
    plot(v_pts(:,2), v_pts(:,1), 'y--');
    plot(v_pts(1,2), v_pts(1,1), 'rx');
    plot(v_pts(end,2), v_pts(end,1), 'gx');
    
    fov_mask = true(size(vessel_patch));
    fov_mask(:,[1 end]) = 0;
    fov_mask([1 end],:) = 0;
    clim = prctile(vessel_patch_equalised(:), [1 99]);
    set(gca, 'CLim', clim);
    cdata = get(i1, 'CData');
    
    loop = true;
    while loop
        selection = questdlg('Add mask?',...
                         'Warning',...
                         'Yes','No','Yes');      

        if strcmpi(selection, 'Yes')
            [ignore_mask] = roipoly;
            fov_mask(ignore_mask') = 0;
            cdata(ignore_mask) = clim(2);
            set(i1, 'CData', cdata);
            set(gca, 'CLim', clim);
        else
            loop = false;
        end
        
    end
      
    save(save_name, 'fov_mask');          

end

close(seg_fig);
display(['Manual segmentation terminated at image ', num2str(i_ve)]);
            