DECOMP_TYPE="{'g2da','h2da'}" NUM_ANGLES=6 WIN_SIZE=3 SIGMA_RANGE="[4 8 16 32 64]" OUTPUT_TYPE="orientation" PREDICTION_TYPE="rf_regression" MASK_DIR="vessel_centre_masks" NUM_TREES=10 IMAGE_TYPE="real" DATA_ROOT="scratch/nailfold/" IMAGE_ROOT="data/set12g" NUM_SAMPLES=50000 MODEL_ROOT="models/vessel" qsub -V -l highmem -hold_ji 258977 -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
DATA_ROOT="scratch/nailfold/models/" MODEL_ROOT="vessel/orientation/rf_regression" MODEL_PATH="259076" MAKE_SAMPLED_MAPS=1 qsub -V -l short -hold_jid 261471 matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
DATA_ROOT="scratch/nailfold/" MODEL_ROOT="models/vessel/orientation/rf_regression" MODEL_PATH="259076" NUM_JOBS=60 IMAGE_ROOT="data/set12g" USE_SAMPLED_MAPS=1 MASK_DIR="fov_masks" qsub -V -t 31-60 -l twoday -hold_jid 261841 matlab_code/trunk/hydra/cuc/predict_image_set.sh
%%
aam_base_dir = 'C:\isbe\nailfold\data\rsa_study\test\aam\';
aam_dirs = dir([aam_base_dir '*c']);
im_tag_len = length('%5Cimages%5C') + 1;
pt_tag_len = length('points%5C') + 1;

for i_dir = 201:length(aam_dirs)
    images_to_move = dir([aam_base_dir aam_dirs(i_dir).name '\%5Cimages%5C*']);
    for i_im = 1:length(images_to_move)
        movefile(...
            [aam_base_dir aam_dirs(i_dir).name '\' images_to_move(i_im).name],...
            [aam_base_dir aam_dirs(i_dir).name '\images\' images_to_move(i_im).name(im_tag_len:end)]);
    end
    
    points_to_move = dir([aam_base_dir aam_dirs(i_dir).name '\points%5C*']);
    for i_im = 1:length(images_to_move)
        movefile(...
            [aam_base_dir aam_dirs(i_dir).name '\' points_to_move(i_im).name],...
            [aam_base_dir aam_dirs(i_dir).name '\points\' points_to_move(i_im).name(pt_tag_len:end)]);
    end
    
end