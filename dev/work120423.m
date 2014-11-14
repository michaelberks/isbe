NUM_TREES=10 PREDICTION_TYPE="orientation" IMAGE_DIR="retinograms/DRIVE/training" ORI_DIR="predictions/orientation/g2d/analytic/orientations_c" SAMPLING_METHOD="generate_vessel_data" WIN_SIZE=3 DO_UBOUND=0 SPLIT_CRITERION="ssq" VAR_CRITERION="ssq" END_CUT_MIN=0 DECOMP_TYPE="dt" NUM_LEVELS=5 qsub -t 1-20 -V -l twoday matlab_code/trunk/hydra/cuc/build_vessel_predictor.sh
qstat
FOREST_JOB="'89706'" FOREST_DIR="models/vessel/orientation" qsub -V -l short -hold_jid 89706 matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
NUM_JOBS=20 FOREST_JOB="'90002'" TEST_IMAGE_DIR="'retinograms/DRIVE/test/images'" MASK_DIR="retinograms/DRIVE/test/foveal_masks" FOREST_DIR="models/vessel/orientation" qsub -l twoday -t 1-20 -V matlab_code/trunk/hydra/cuc/classify_image_set.sh
qstat
%%

mkdir([asymmetryroot('shared') 'data\mammograms\mass_roi\orientations']);


for ii = 1:146
    %load retinogram and merge RGB channels
    mammo = u_load([asymmetryroot('shared') 'data\mammograms\mass_roi\images\bg' zerostr(ii,5) '.mat']);

    %Compute repsonses for mono and save
    [response_map ori_map scale_map] = gaussian_2nd_derivative_line(mammo, [1.1 1.9 3.2]);
    save_uint8([asymmetryroot('shared') 'data\mammograms\mass_roi\orientations\ori' zerostr(ii,5) '.mat'], complex(cos(2*ori_map), sin(2*ori_map)));
    
    mask = u_load([asymmetryroot('shared') 'data\mammograms\mass_roi\foveal_masks\bg' zerostr(ii,5) '_mask.mat']);
    mask([1 end], 1:end) = 0;
    mask(1:end, [1 end]) = 0;
    save([asymmetryroot('shared') 'data\mammograms\mass_roi\foveal_masks\bg' zerostr(ii,5) '_mask.mat'], 'mask');
end
%%
for ii = 121:146
    %load retinogram and merge RGB channels
    mammo = u_load([asymmetryroot('shared') 'data\mammograms\mass_roi\images\bg' zerostr(ii,5) '.mat']);
    mask = u_load([asymmetryroot('shared') 'data\mammograms\mass_roi\foveal_masks\bg' zerostr(ii,5) '_mask.mat']);
    
    figure; 
    subplot(1,2,1); imgray(mammo);
    subplot(1,2,2); imgray(mask);
end