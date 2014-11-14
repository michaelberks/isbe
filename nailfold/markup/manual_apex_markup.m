function manual_apex_markup(nailfold_dir, vessel_dir, start_num, file_type)

if nargin < 3
    start_num = 1;
end
if nargin < 4
    file_type = '.txt';
end

nailfold_list = dir([vessel_dir, '*', file_type]);

seg_fig = figure;
set(seg_fig, 'WindowButtonDownFcn', @get_pts);

for ii = start_num:length(nailfold_list)
    image_name = [nailfold_dir, nailfold_list(ii).name(1:end-11) '.bmp'];
    
    if strcmp(file_type, '.mat')
        nailfold = u_load(image_name);
    else
        nailfold = imread(image_name);
    end
    
    %Load in the annotated vessels
    v = read_vessels_from([vessel_dir nailfold_list(ii).name]);
    num_vessels = length(v);
    [rows cols] = size(nailfold);
    
    vessel_pts = zeros(num_vessels, 3);
    
    for jj = 1:num_vessels       
        %Discard duplicate points
        keep = [true; any(diff(v{jj}),2)];
        vessel = [v{jj}(keep,1) v{jj}(keep,2)];

        %Sample points evenly along the vessel
        dists = cumsum([0; sqrt(sum(diff(vessel).^2, 2))]);
        vessel = interp1(dists, vessel, linspace(0, dists(end), floor(dists(end))), 'spline');
        
        %Get top of vessel as min y-value
        [y_apex min_idx] = min(vessel(:,2));
        x_apex = vessel(min_idx,1);
        
        %Workout ROI in image to extract
        sr = max(1, floor(y_apex - 25));
        er = min(rows, ceil(y_apex + 75));
        sc = max(1, floor(x_apex - 100));
        ec = min(cols, floor(x_apex + 100));
        image_patch = nailfold(sr:er, sc:ec);
        
        %Correct vessel coords to ROI frame
        vessel(:,1) = vessel(:,1) - sc;
        vessel(:,2) = vessel(:,2) - sr;
        
        %Display image patch and plot vessel
        figure(seg_fig);
        imgray(image_patch);
        plot(vessel(:,1), vessel(:,2), 'r.', 'markersize', 2);
        set(seg_fig, 'Name', ['Vessel ', num2str(jj), ' of ', num2str(num_vessels)]);

        %Select start of top section
        [xi,yi,~] = impixel;
        dists = (vessel(:,1)-xi).^2 + (vessel(:,2)-yi).^2;
        [~, v_start] = min(dists);
        plot(vessel(v_start,1), vessel(v_start,2), 'go');
        
        %Select apex
        [xi,yi,~] = impixel;
        dists = (vessel(:,1)-xi).^2 + (vessel(:,2)-yi).^2;
        [~, v_apex] = min(dists);
        plot(vessel(v_apex,1), vessel(v_apex,2), 'go');
        
        %Select end of top
        [xi,yi,~] = impixel;
        dists = (vessel(:,1)-xi).^2 + (vessel(:,2)-yi).^2;
        [~, v_end] = min(dists);
        plot(vessel(v_end,1), vessel(v_end,2), 'go');
        
        %Save marked points
        vessel_pts(jj,:) = [v_start v_apex v_end];
        save([vessel_dir nailfold_list(ii).name(1:end-4) '_v_pts.mat'], 'vessel_pts');
        delete(get(seg_fig, 'children'));
        
    end
    clear nailfold;
    
    
    set(0,'DefaultFigureWindowStyle','normal');
    selection = questdlg('Continue to next image?',...
                     'Warning',...
                     'Yes','No','Yes');
    set(0,'DefaultFigureWindowStyle','docked');             
    if strcmpi(selection, 'No')
        break
    end
    
end

close(seg_fig);
display(['Manual segmentation terminated at image ', num2str(ii)]);
    
            