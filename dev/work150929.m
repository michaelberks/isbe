load('C:\isbe\nailfold\data\rsa_study\data_lists\miccai_lists.mat', 'miccai_selection');
load('C:\isbe\nailfold\data\rsa_study\data_lists\image_id_data.mat');
im_names = sort(image_id_data.im_names(miccai_selection.validation));

for i_im = 1:length(im_names)
    im_name = im_names{i_im};
    display(['Processing image ' num2str(i_im) ': ' im_name]);
    
    root_dir = 'C:\isbe\nailfold\data\rsa_study\master_set\';

    %vessel_predictions = load([root_dir 'detected_capillaries_cxx\' im_name '_vessels_v_pred.txt']);
    cxx_candidates = load([root_dir 'detected_capillaries_cxx\' im_name '_apex_candidates.txt']);
    candidate_xy = cxx_candidates(:,1:2)+1;
    candidate_scores = cxx_candidates(:,6);
    candidate_oris = complex(cxx_candidates(:,4), cxx_candidates(:,5));
    candidate_widths = cxx_candidates(:,3);
    candidate_rescores = cxx_candidates(:,7);
    candidate_displacements = cxx_candidates(:,8);

    candidate_class_probs = cxx_candidates(:,10);
    candidate_class = cxx_candidates(:,9);
    
    save_path = ['C:\isbe\nailfold\data\rsa_study\master_set\detected_capillaries_cxx\mat\' im_name '_caps.mat'];
    
    
    if ~isempty(candidate_class) && candidate_class(1) >= 0
        status = 0;
%         load(save_path,...
%             'selected_distal', 'selected_non_distal', 'merged_with');
        
        load(save_path,...
            'merged_with');
        
        selected_distal = candidate_class == 1;
        selected_non_distal = candidate_class == 2;
% 
%         [merged_with] = post_merge_candidates(candidate_xy, candidate_rescores, candidate_widths, ...
%              vessel_predictions(:,:,1), 60, 0.5, 20, 0);
%         
% 
%         for i_can = 1:length(merged_with)
% 
%             if merged_with(i_can)
% 
%                 if selected_distal(i_can)
% 
%                     if selected_distal(merged_with(i_can))
%                         %Two distals merged together, discard the one with
%                         %lower score
%                         selected_distal(i_can) = 0;
% 
%                     elseif selected_non_distal(merged_with(i_can))
%                         %Distal, merged with a non-distal. The
%                         %non-distal has higher score, so must have been
%                         %reject on its location - do nothing
% 
%                     else
%                         %Distal, merged with something not selected as
%                         %a non-distal, but has a higher vessel score
%                         %than it. Swap the distal status
%                         selected_distal(i_can) = 0;
%                         selected_distal(merged_with(i_can)) = 1;
%                     end
%                 elseif selected_non_distal(i_can)
% 
%                     if selected_distal(merged_with(i_can))
%                         %Non-distal, merged with a distal. Discard the
%                         %non-distal altogether
%                         selected_non_distal(i_can) = 0;
% 
%                     elseif selected_non_distal(merged_with(i_can))
%                         %2 non-distals merged together, keep the one
%                         %with the higher score
%                         selected_non_distal(i_can) = 0;
% 
%                     else
%                         %Non-istal, merged with something not selected as
%                         %a non-distal, but has a higher vessel score
%                         %than it. Swap the distal status
%                         selected_non_distal(i_can) = 0;
%                         selected_non_distal(merged_with(i_can)) = 1;
%                     end
%                 end
%             end
%         end


        fell_at_the_last = false(size(candidate_rescores));
        %Finally, we discard any vessels that lie on top of one another
        if any(selected_distal)
            %Which candidates belong to the distal row?
            [distal_idx] = select_distal_candidates(...
                candidate_xy(selected_distal,:), candidate_rescores(selected_distal,:),...
                10, 1.309, 0, 0);

            %Work out which ones we've discarded (with indices relative to the
            %original list)
            distal_idx_tf = true(sum(selected_distal),1);
            distal_idx_tf(distal_idx) = 0;
            fell_at_the_last(selected_distal) = distal_idx_tf;

            %
            selected_distal(fell_at_the_last) = 0;
            selected_non_distal(fell_at_the_last & candidate_displacements > 0) = 1;
        end
    else
        %Minimum number of candidates not met to estimate displacements,
        %apply strong threshold
        selected_distal = candidate_rescores > 0.8;
        selected_non_distal = false(size(candidate_rescores));
        fell_at_the_last = [];
        merged_with = [];
        status = -1;
    end

    save(save_path,...
            'selected_distal', 'selected_non_distal', 'merged_with', 'status');    
    
end