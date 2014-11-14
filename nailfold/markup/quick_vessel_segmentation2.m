function quick_vessel_segmentation2(vessel_dir, vessel_selection)



vessel_list = dir([vessel_dir, '*vessel.mat']);


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
    vessel_name = [vessel_dir, vessel_list(i_ve).name];
    
    vessel_struc = u_load(vessel_name);
    vessel_patch = vessel_struc.vessel_patch;
    v_pts = vessel_struc.v_pts;
    
    
    figure(seg_fig);
    set(seg_fig, 'Name', ['Vessel ', num2str(counter), ' of ', num2str(num_vessels),...
        ': ', vessel_list(i_ve).name]);
    
    subplot(1,2,2);
    imgray(vessel_patch);
    plot(v_pts(:,1), v_pts(:,2));
  
    selection = questdlg('How good is it really?',...
                     'Warning',...
                     'Ok','Great','Ok');
                 
    switch selection
        case 'Great'            
            vessel_struc.quality = 3;
            save(vessel_name, 'vessel_struc');
    end          
    
    selection = questdlg('Continue to next image?',...
                     'Warning',...
                     'Yes','No','Yes');           
    if strcmpi(selection, 'No')
        break
    end
    delete(get(seg_fig, 'children'));
end

close(seg_fig);
display(['Manual segmentation terminated at image ', num2str(i_ve)]);
            