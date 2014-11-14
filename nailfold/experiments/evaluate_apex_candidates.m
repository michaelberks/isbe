function [candidate_label gt_hits candidate_match] = evaluate_apex_candidates(gt_xy, candidate_xy, tol, ...
    vessel_prob, connectivity_tol, connectivity_dist, do_plot)


if ~exist('connectivity_tol', 'var') || isempty(connectivity_tol)
    connectivity_tol = 0.5;
end
if ~exist('connectivity_dist', 'var') || isempty(connectivity_dist)
    connectivity_dist = 100;
end
if ~exist('tol', 'var') || isempty(tol)
    tol = 15;
end
if ~exist('vessel_prob', 'var')
    vessel_prob = [];
end
if ~exist('do_plot', 'var')|| isempty(do_plot)
    do_plot = 0;
end
patch_sz2 = ceil(connectivity_dist / 2) + 1;
patch_sz = patch_sz2*2 + 1;
connectivity_dist = connectivity_dist^2;
plot_num = 1;
n_connect_pts = 20;

%First work out, is candidate true or false?

num_cans = size(candidate_xy, 1);
num_gt = size(gt_xy, 1);

if numel(tol) == 1
    tol = ones(num_gt,1)*tol;
end
tol2 = tol.^2;

candidate_label = zeros(num_cans, 1);
gt_hits = zeros(num_gt, 1);
candidate_match = zeros(num_cans, 1);

if ~num_gt
    return;
end

for i_can = 1:num_cans
    
    [min_dist gt_idx] = min( sum( bsxfun(@minus, gt_xy, candidate_xy(i_can,:)).^2, 2) );
    candidate_label(i_can) = min_dist < tol2(gt_idx);
    
    if ~candidate_label(i_can) && ~isempty(vessel_prob) && (min_dist < connectivity_dist)
        
        %Get the two centres and their midpoint
        centre1 = candidate_xy(i_can,:);
        centre2 = gt_xy(gt_idx,:);
        centre_midpoint = (centre1 + centre2) / 2;

        %Sample a patch about the midpoint
        vessel_prob_patch = sample_window(vessel_prob, patch_sz, ...
            round(centre_midpoint(2)), round(centre_midpoint(1)));

        %Make the centres xy coords relative to the patch frame
        centre1 = centre1 - round(centre_midpoint) + patch_sz2;
        centre2 = centre2 - round(centre_midpoint) + patch_sz2;    

        %Work connectivity
        for i_con = fliplr(linspace(0, 1, n_connect_pts));

            con_mask = bwselect(vessel_prob_patch > i_con, centre1(1), centre1(2));
            if con_mask(round(centre2(2)), round(centre2(1)))                            
                break;
            end
        end
        connectedness = i_con;
        
       candidate_label(i_can) = connectedness > connectivity_tol;
       
       if do_plot
            %Display
            if plot_num == 1
                figure;
            end
            subplot(3,4,plot_num); 
            imgray(vessel_prob_patch); 

            xlabel(num2str(connectedness));

            plot(centre1(:,1), centre1(:,2), 'x', 'markersize', 20);
            plot(centre2(:,1), centre2(:,2), 'x', 'markersize', 20);
            
            %If we discarded either candidate circle the one we kept
            if candidate_label(i_can)
                plot(centre1(:,1), centre1(:,2), 'ro', 'markersize', 20);
            end

            if plot_num == 12
                plot_num = 1;
            else
                plot_num = plot_num + 1;
            end
        end
           
    end
        
        
    if candidate_label(i_can)
        gt_hits(gt_idx) = i_can;
        candidate_match(i_can) = gt_idx;
    end
    
end


% num_gt = size(gt_xy, 1);
% tp_count = cumsum(candidate_label);
% 
% fp_count = (1:num_cans)' - tp_count;
% 
% sensitivity = tp_count / num_gt;
% 
% figure; plot(fp_count, sensitivity);
    
    



