%First of all get the dat asomewhere sensible and make some vessel masks

nf_files = dir('C:\isbe\nailfold\images\anonymous_oct\annotations_qt\*.txt');
num_nf = length(nf_files);
%%

mkdir('C:\isbe\nailfold\data\training\images');
mkdir('C:\isbe\nailfold\data\training\fov_masks');
mkdir('C:\isbe\nailfold\data\training\vessel_masks');
mkdir('C:\isbe\nailfold\data\test\images');
mkdir('C:\isbe\nailfold\data\test\fov_masks');
mkdir('C:\isbe\nailfold\data\test\vessel_masks');
mkdir('C:\isbe\nailfold\data\training2\images');
mkdir('C:\isbe\nailfold\data\training2\fov_masks');
mkdir('C:\isbe\nailfold\data\training2\vessel_masks');
mkdir('C:\isbe\nailfold\data\test2\images');
mkdir('C:\isbe\nailfold\data\test2\fov_masks');
mkdir('C:\isbe\nailfold\data\test2\vessel_masks');
%%
do_plot = 0;
do_half = 0;
for nn = 1:12
    
    %Load in nailfold
    nailfold = imread(['C:\isbe\nailfold\images\anonymous_oct\bmp\' nf_files(nn).name(1:end-11) '.bmp']);
    nailfold = nailfold(:,:,1);
    if do_half
        nailfold = imresize(nailfold, 0.5, 'bilinear');
    end
    [rows cols] = size(nailfold);
    
    fov_mask = imerode(nailfold > 10 & nailfold < 250, strel('disk', 10));
    vessel_mask = false(rows, cols);
  
    %Load in the annotated vessels
    v = read_vessels_from(['C:\isbe\nailfold\images\anonymous_oct\annotations_qt\' nf_files(nn).name]);
    num_vessels = length(v);
    
    %Add each vessel to the vessel mask
    for ii = 1:num_vessels
       
        if do_half
            ri = round(v{ii}(:,2)/2);
            ci = round(v{ii}(:,1)/2);
        else
            ri = round(v{ii}(:,2));
            ci = round(v{ii}(:,1));
        end
            
        v_idx = sub2ind([rows cols], ri, ci);
        vessel_mask(v_idx) = 1;
    end
    
    %Dilate the vessel mask to widen the vessels
    vessel_mask = imdilate(vessel_mask, strel('disk', 1));
    
    if do_plot
        figure;
        subplot(1,2,1); imgray(nailfold);
        subplot(1,2,2); imgray(vessel_mask);
    end
    
    if nn < 8
        type = 'training';
    else
        type = 'test';
    end
    if do_half
        root = ['C:\isbe\nailfold\data\' type '2'];
    else
        root = ['C:\isbe\nailfold\data\' type];
    end
    save([root '\images\' zerostr(nn,2) '_' type '.mat'], 'nailfold');
    save([root '\fov_masks\' zerostr(nn,2) '_' type '_fov_mask.mat'], 'fov_mask');
    save([root '\vessel_masks\' zerostr(nn,2) '_' type '_vessel_mask.mat'], 'vessel_mask');
    
end
display(num_pts);
%%
root = 'C:\isbe\nailfold\data\training';
for nn = 1:7
    nailfold = u_load([root '\images\' zerostr(nn,2) '_training.mat']);
    fov_mask = u_load([root '\fov_masks\' zerostr(nn,2) '_training_fov_mask.mat']);
    vessel_mask = u_load([root '\vessel_masks\' zerostr(nn,2) '_training_vessel_mask.mat']);
    
    c_min = min(nailfold(fov_mask));
    c_max = max(nailfold(fov_mask));
    
    figure;
    a1 = subplot(2,1,1); imgray(nailfold); caxis([c_min c_max]);
    a2 = subplot(2,1,2); imgray(vessel_mask); 
    linkaxes([a1 a2]);
    %plot(vx, vy, 'r.', 'markersize', 2);
end
%%
root = 'C:\isbe\nailfold\data\test';
for nn = 8:12
    vessel_prob = load_uint8(['C:\isbe\nailfold\data\test\images\results\107629\' zerostr(nn,2) '_test_class.mat']);
    nailfold = u_load([root '\images\' zerostr(nn,2) '_test.mat']);
    fov_mask = u_load([root '\fov_masks\' zerostr(nn,2) '_test_fov_mask.mat']);
    vessel_mask = u_load([root '\vessel_masks\' zerostr(nn,2) '_test_vessel_mask.mat']);
    
    vessel_prob(vessel_mask) = 0;
    
    c_min = min(nailfold(fov_mask));
    c_max = max(nailfold(fov_mask));
    
    [vy vx] = find(vessel_mask);
    figure;
    subplot(2,1,1); imgray(nailfold); caxis([c_min c_max]);
    subplot(2,1,2); imgray(vessel_prob); 
    %plot(vx, vy, 'r.', 'markersize', 2);
end
%%-------------------------------------------------------------------------
%Round 2:
%%
contour_dir = 'C:\isbe\nailfold\data\rsa_study\vessel_contours\normal\';
vessel_dir = 'C:\isbe\nailfold\data\rsa_study\apexes\normal\';
vessel_mask_dir = 'C:\isbe\nailfold\data\rsa_study\training\vessel_masks\';
vessel_image_dir = 'C:\isbe\nailfold\data\rsa_study\training\images\';
total_v_pts = 0;
for i_v = 1:length(v_files)
    load([vessel_mask_dir v_files(i_v).name(1:end-12) '_v_mask.mat'], 'vessel_mask');
    total_v_pts = total_v_pts + sum(vessel_mask(:));
end
%%
DECOMP_TYPE="{'g2da','h2da'}" NUM_ANGLES=6 WIN_SIZE=3 SIGMA_RANGE="[4 8 16]" OUTPUT_TYPE="detection" PREDICTION_TYPE="rf_classification" NUM_TREES=2 IMAGE_TYPE="real" DATA_ROOT="scratch/nailfold/data" IMAGE_ROOT="training" NUM_SAMPLES=5000 MODEL_ROOT="models/vessel" qsub -V -l twoday -t 1-20 matlab_code/trunk/hydra/cuc/build_predictor.sh
MODEL_ROOT="models/vessel/detection/rf_classification" MODEL_PATH="'88126'" qsub -V -hold_jid 88126 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/vessel/detection/rf_classification" MODEL_PATH="'88126'" NUM_JOBS=10 IMAGE_ROOT="retinograms/DRIVE/test" MASK_DIR="fov_masks" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 88127 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
qstat
