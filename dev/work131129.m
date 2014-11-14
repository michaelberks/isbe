DECOMP_TYPE="{'g2dia','h2dia'}" NUM_ANGLES=6 WIN_SIZE=3 SIGMA_RANGE="[2 5]" OUTPUT_TYPE="detection" PREDICTION_TYPE="rf_classification" MASK_DIR="vessel_centre_masks" NUM_TREES=1 IMAGE_TYPE="real" DATA_ROOT="scratch/nailfold/" IMAGE_ROOT="data/rsa_study/set12g_half" NUM_SAMPLES=50000 MODEL_ROOT="models/vessel" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="{'g2dia','h2dia'}" NUM_ANGLES=6 WIN_SIZE=3 SIGMA_RANGE="[4 5]" OUTPUT_TYPE="detection" PREDICTION_TYPE="rf_classification" MASK_DIR="vessel_centre_masks" NUM_TREES=1 IMAGE_TYPE="real" DATA_ROOT="scratch/nailfold/" IMAGE_ROOT="data/rsa_study/set12g" NUM_SAMPLES=5000 MODEL_ROOT="models/vessel" qsub -V -l short -t 1 matlab_code/trunk/hydra/cuc/build_predictor.sh
create_image_lists
generate_training_data
compute_filter_responses

DECOMP_TYPE="{'g2dia','h2dia'}" NUM_ANGLES=6 WIN_SIZE=3 SIGMA_RANGE="[2 5]" OUTPUT_TYPE="orientation" PREDICTION_TYPE="rf_regression" MASK_DIR="vessel_centre_masks" NUM_TREES=10 IMAGE_TYPE="real" DATA_ROOT="scratch/nailfold/" IMAGE_ROOT="data/rsa_study/set12g_half" NUM_SAMPLES=50000 MODEL_ROOT="models/vessel" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="{'g2dia','h2dia'}" NUM_ANGLES=6 WIN_SIZE=3 SIGMA_RANGE="[2 5]" OUTPUT_TYPE="width" PREDICTION_TYPE="rf_regression" MAX_N_IMAGES=450 SHIFT_IMAGES=450 MASK_DIR="vessel_centre_masks" NUM_TREES=10 IMAGE_TYPE="real" DATA_ROOT="scratch/nailfold/" IMAGE_ROOT="data/rsa_study/set12g_half" NUM_SAMPLES=50000 MODEL_ROOT="models/vessel" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh

DATA_ROOT="scratch/nailfold/models/" MODEL_ROOT="vessel/detection/rf_classification" MODEL_PATH="294370" MAKE_SAMPLED_MAPS=1 qsub -V -l short -hold_jid 294370 matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
DATA_ROOT="scratch/nailfold/models/" MODEL_ROOT="vessel/orientation/rf_regression" MODEL_PATH="294378" MAKE_SAMPLED_MAPS=1 qsub -V -l short -hold_jid 294378 matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
DATA_ROOT="scratch/nailfold/models/" MODEL_ROOT="vessel/width/rf_regression" MODEL_PATH="294379" MAKE_SAMPLED_MAPS=1 qsub -V -l short -hold_jid 294379 matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh

DATA_ROOT="scratch/nailfold/" MODEL_ROOT="models/vessel/detection/rf_classification" MODEL_PATH="294672" NUM_JOBS=196 IMAGE_ROOT="data/rsa_study/test_half" USE_SAMPLED_MAPS=0 OVERWRITE=0 qsub -V -t 1-10 matlab_code/trunk/hydra/cuc/predict_image_set.sh

%%
hole_gs = [];
can_labels = [];
load('C:\isbe\nailfold\data\rsa_study\test\results\rf_offsets_20131121T181011.mat');
 
for i_im = 1:length(candidates_list)
    
    display(['processing image ' num2str(i_im)]);
    
    im_name = candidates_list(i_im).name(1:6);
    
    load(['C:\isbe\nailfold\data\rsa_study\test\apex_maps\frog\local_maxima2\' candidates_list(i_im).name]);
    
    num_candidates = length(candidate_scores);
    
    for i_can = 1:num_candidates
        
        load(['C:\isbe\nailfold\data\rsa_study\test\aam\' im_name '\apex' zerostr(i_can,4) '.mat']);
        if i_can == 1
            
            nailfold = u_load(['C:\isbe\nailfold\data\rsa_study\test\images\' im_name '.mat']);
            
        end
        
        vessel_patch = nailfold(apex_candidate.sr:apex_candidate.er, apex_candidate.sc:apex_candidate.ec);
        
        hole_mask = poly2mask(apex_candidate.vessel_xy(:,1), apex_candidate.vessel_xy(:,2), size(vessel_patch,1), size(vessel_patch,2));
        
        hole_g = mean(vessel_patch(hole_mask));
        
        hole_gs = [hole_gs; hole_g];
        
    end
    
    can_labels = [can_labels; detections{i_im,3} > 0];
end
%%
for i_im = 1:length(candidates_list)
    
    display(['processing image ' num2str(i_im)]);
    
    im_name = candidates_list(i_im).name(1:6);
    
    load(['C:\isbe\nailfold\data\rsa_study\test\apex_maps\frog\local_maxima2\' candidates_list(i_im).name]);
    
    num_candidates = length(candidate_scores);
    
    if any(kept)
        nailfold = u_load(['C:\isbe\nailfold\data\rsa_study\test\images\' im_name '.mat']);

        for i_can = 1:num_candidates    
            if kept(i_can)
                load(['C:\isbe\nailfold\data\rsa_study\test\aam\' im_name '\apex' zerostr(i_can,4) '.mat']);
                vessel_patch = nailfold(apex_candidate.sr:apex_candidate.er, apex_candidate.sc:apex_candidate.ec);

                if detections{i_im,3}(i_can)
                    write_im_from_colormap(vessel_patch, ['C:\isbe\nailfold\data\rsa_study\test\results\misclassified_apexes\tps\' im_name '_apex' zerostr(i_can,4) '.png']);
                else
                    write_im_from_colormap(vessel_patch, ['C:\isbe\nailfold\data\rsa_study\test\results\misclassified_apexes\fps\' im_name '_apex' zerostr(i_can,4) '.png']);
                end
            end    
        end
    end

end
%%
mkdir C:\isbe\nailfold\data\rsa_study\test_half\images;
for i_im = 1:length(candidates_list)
    im_name = candidates_list(i_im).name(1:6);
    nailfold = u_load(['C:\isbe\nailfold\data\rsa_study\test\images\' im_name '.mat']);
    nailfold = imresize(nailfold, 0.5);
    save(['C:\isbe\nailfold\data\rsa_study\test_half\images\' im_name '.mat'], 'nailfold');
end
        