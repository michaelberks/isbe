function quick_vessel_segmentation(vessel_dir, start_num)

if nargin < 2
    start_num = 1;
end

vessel_list = dir([vessel_dir, '*vessel.mat']);
screen_size = get(0,'ScreenSize');
seg_fig = figure(...
    'windowstyle', 'normal',...
    'position', screen_size);
for i_ve = start_num:length(vessel_list)
    vessel_name = [vessel_dir, vessel_list(i_ve).name];
    
    vessel_struc = load(vessel_name);
    vessel_patch = vessel_struc.vessel_patch;
    v_pts = vessel_struc.v_pts;
    
    
    figure(seg_fig);
    set(seg_fig, 'Name', ['Vessel ', num2str(i_ve), ' of ', num2str(length(vessel_list)),...
        ': ', vessel_list(i_ve).name]);
    
    subplot(1,2,2);
    imgray(vessel_patch);
    plot(v_pts(:,1), v_pts(:,2));
  
    selection = questdlg('Continue to next image?',...
                     'Warning',...
                     'Reject','Skip','Process','Process');
                 
    switch selection
        case 'Reject'            
            vessel_struc.quality = 0;

        case 'Skip'
            vessel_struc.quality = 1;

        case 'Process'
            vessel_struc.quality = 2;
    end
    save(vessel_name, 'vessel_struc');
    
    
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
            