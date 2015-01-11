num_cells = size(cell_positions,1);
colors = lines(num_cells);
figure;
for i_p = 1:9
    subplot(3,3, i_p); hold on;
    for i_c = 1:num_cells
        plot(cell_positions(i_c,1,i_p), cell_positions(i_c,2,i_p), '.', 'markerfacecolor', colors(i_c,:));
        text(cell_positions(i_c,1,i_p), cell_positions(i_c,2,i_p), num2str(i_c));
        
    end
end
%%
delete('C:\isbe\nailfold\synthesis\showcase\ncm_outline\input\logo.gif');
map = gray(256);
for i_f = 1:100
    frame = imread(['C:\isbe\nailfold\synthesis\showcase\ncm_outline\input\frame_' zerostr(i_f, 4) '.png']);

    if i_f == 1
        imwrite(frame, map, 'C:\isbe\nailfold\synthesis\showcase\ncm_outline\input\logo.gif', 'gif', 'WriteMode', 'overwrite', 'Loop', Inf, 'DelayTime', 0.01);
    else
        imwrite(frame, map, 'C:\isbe\nailfold\synthesis\showcase\ncm_outline\input\logo.gif', 'gif', 'WriteMode', 'append', 'DelayTime', 0.01);
    end
end
%%
frame = imread(['C:\isbe\nailfold\synthesis\showcase\ncm_outline\input\frame_' zerostr(1, 4) '.png']);
imsz = [588 620];

sr = round(linspace(1,15,15));
er = round(linspace(174,imsz(1),15));
sc = round(linspace(1,215,15));
ec = round(linspace(354,imsz(2),15));

figure;
for i_f = 1:15
frame_i = imresize(frame(sr(i_f):er(16-i_f), sc(i_f):ec(16-i_f)),[200 200]);
    subplot(3,5,i_f); imgray(frame_i);
end
%%
imsz = [248 320];

sr = [ones(1,25) round(linspace(1,14,50)) 25*ones(1,25)];
er = [58*ones(1,25) round(linspace(58,imsz(1),50)) imsz(1)*ones(1,25)];
sc = [ones(1,25) round(linspace(1,180,50)) 180*ones(1,25)];
ec = [224*ones(1,25) round(linspace(224,imsz(2),50)) imsz(2)*ones(1,25)];

output_dir = 'C:\isbe\nailfold\synthesis\showcase\bms\logo4\';
logoname = fullfile(output_dir, 'logo_pink_zoom.gif');
map = pink(256);

for i_f = 1:200
    frame = imread(fullfile(output_dir, ['frame_' zerostr(i_f, 4) '.png']));
    frame = double(frame);
    frame = 255*(frame - min(frame(:))) / (max(frame(:))-min(frame(:)));
    frame = uint8(frame);

    if i_f < 101
        frame_i = imresize(frame(sr(101-i_f):er(i_f), sc(101-i_f):ec(i_f)),[186 240]);
    else
        frame_i = imresize(frame(sr(i_f-100):er(201-i_f), sc(i_f-100):ec(201-i_f)),[186 240]);
    end

    if i_f == 1
        imwrite(frame_i, map, logoname, 'gif', 'WriteMode', 'overwrite', 'Loop', Inf, 'DelayTime', 0.01);
    else
        imwrite(frame_i, map, logoname, 'gif', 'WriteMode', 'append', 'DelayTime', 0.01);
    end
end
%%
bg = imread('C:\isbe\nailfold\synthesis\showcase\bms\bg4.png');

sr1 = 360;
sc1 = 830;
sr2 = 352;
sc2 = 1513;

bg_roi1 = bg(sr1+(1:320), sc1+(1:320));
bg_roi2 = bg(sr2+(1:320), sc2+(1:320));

bms_outline_video('P:\isbe\nailfold\bms_outline2.png', 'C:\isbe\nailfold\synthesis\showcase\bms\logo6', bg_roi1, 0.2, 42, 2.0);
bms_outline_video('P:\isbe\nailfold\bms_outline5.png', 'C:\isbe\nailfold\synthesis\showcase\bms\logo7', bg_roi2, 0.2, 42, 1.0);

logo1 = imread('C:\isbe\nailfold\synthesis\showcase\bms\logo6\logo.png');
logo2 = imread('C:\isbe\nailfold\synthesis\showcase\bms\logo7\logo.png');

bg(sr1+(1:size(logo1,1)), sc1+(1:320)) = logo1;
bg(sr2+(1:size(logo2,1)), sc2+(1:320)) = logo2;

figure; imgray(bg);
%%
bg = imread('C:\isbe\nailfold\synthesis\showcase\bms\bg3.png');
sr1 = 260;
sc1 = 1040;
bg_roi = bg(sr1+(1:320), sc1+(1:320));
bms_outline_video('P:\isbe\nailfold\bms_outline2.png', 'C:\isbe\nailfold\synthesis\showcase\bms\logo4', bg_roi, 0.2, 42, 2.0);

logo = imread('C:\isbe\nailfold\synthesis\showcase\bms\logo4\logo.png');
bg(sr1+(1:size(logo,1)), sc1+(1:320)) = logo;
figure; imgray(bg); colormap(pink(256));
%%
output_dir = 'C:\isbe\nailfold\synthesis\showcase\bms\logo7\';
logoname = fullfile(output_dir, 'logo_pink.gif');
map = pink(256);

for i_f = 1:200
    frame = imread(fullfile(output_dir, ['frame_' zerostr(i_f, 4) '.png']));
    frame = double(frame);
    frame = 255*(frame - min(frame(:))) / (max(frame(:))-min(frame(:)));
    frame = uint8(frame);

    if i_f == 1
        imwrite(frame, map, logoname, 'gif', 'WriteMode', 'overwrite', 'Loop', Inf, 'DelayTime', 0.01);
    else
        imwrite(frame, map, logoname, 'gif', 'WriteMode', 'append', 'DelayTime', 0.01);
    end

end
%%
DECOMP_TYPE="{'dt'}" WIN_SIZE=3 OUTPUT_TYPE="orientation" PREDICTION_TYPE="rf_regression" FG_DIR="vessel_masks" NUM_TREES=10 IMAGE_TYPE="real" DATA_ROOT="scratch/nailfold/" IMAGE_ROOT="data/rsa_study/set12g_half" NUM_SAMPLES=50000 MODEL_ROOT="models/vessel" qsub -V -l twoday -t 1 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="{'dt'}" WIN_SIZE=3 OUTPUT_TYPE="detection" PREDICTION_TYPE="rf_classification" FG_DIR="vessel_masks" NUM_TREES=10 IMAGE_TYPE="real" DATA_ROOT="scratch/nailfold/" IMAGE_ROOT="data/rsa_study/set12g_half" NUM_SAMPLES=50000 MODEL_ROOT="models/vessel" qsub -V -l twoday matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="{'dt'}" WIN_SIZE=3 OUTPUT_TYPE="width" PREDICTION_TYPE="rf_regression" FG_DIR="vessel_masks" NUM_TREES=10 IMAGE_TYPE="real" DATA_ROOT="scratch/nailfold/" IMAGE_ROOT="data/rsa_study/set12g_half" NUM_SAMPLES=50000 MODEL_ROOT="models/vessel" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
qstat
DATA_ROOT="scratch/nailfold/models/" MODEL_ROOT="vessel/orientation/rf_regression" MODEL_PATH="647786" MAKE_SAMPLED_MAPS=1 qsub -V -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
DATA_ROOT="scratch/nailfold/models/" MODEL_ROOT="detection/rf_classification" MODEL_PATH="647785" MAKE_SAMPLED_MAPS=1 qsub -V -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
DATA_ROOT="scratch/nailfold/models/" MODEL_ROOT="vessel/width/rf_regression" MODEL_PATH="647780" MAKE_SAMPLED_MAPS=1 qsub -V -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh

DATA_ROOT="scratch/nailfold/" MODEL_ROOT="models/vessel/width/rf_regression" MODEL_PATH="647780" NUM_JOBS=50 IMAGE_ROOT="data/rsa_study/master_set" USE_SAMPLED_MAPS=0 OVERWRITE=0 qsub -l twoday -V -t 1-50 matlab_code/trunk/hydra/cuc/predict_image_set.sh
DATA_ROOT="scratch/nailfold/" MODEL_ROOT="models/vessel/detection/rf_classification" MODEL_PATH="647785" NUM_JOBS=50 IMAGE_ROOT="data/rsa_study/master_set" USE_SAMPLED_MAPS=0 OVERWRITE=0 qsub -l twoday -V -t 1-50 matlab_code/trunk/hydra/cuc/predict_image_set.sh
DATA_ROOT="scratch/nailfold/" MODEL_ROOT="models/vessel/orientation/rf_regression" MODEL_PATH="647786" NUM_JOBS=50 IMAGE_ROOT="data/rsa_study/master_set" USE_SAMPLED_MAPS=0 OVERWRITE=0 qsub -l twoday -V -t 2-50 matlab_code/trunk/hydra/cuc/predict_image_set.sh

%%
fn = fieldnames(a);

for i_f = 1:length(fn);
    cmd = [fn{i_f} ' = ['];
    
    for i_d = 1:length(a.(fn{i_f}))
        cmd = [cmd num2str(a.(fn{i_f})(i_d)) '; ']; %#ok
    end
    display([cmd(1:end-1) '];']);
end
%%
load('C:\isbe\nailfold\data\rsa_study\data_lists\miccai_lists.mat');
run_ncm_batch('rsa_study/master_set', 1,... % the user's input
    'im_idx',               miccai_selection.test,...
    'maxima_name',          'miccai_maxima',...
    'prob_dir',             'rf_classification/647785/',... 296655
    'ori_dir',              'rf_regression/647786/',...
    'width_dir',            'rf_regression/647780/',...
    'fov_mask_dir',         'fov_masks/',...
    'centre_dir',           'vessel_centres/dt/');
%%
DATA_ROOT="scratch/nailfold/" MODEL_ROOT="models/apex" MODEL_PATH="set12g_half_296655" NUM_JOBS=50 IMAGE_ROOT="data/rsa_study/master_set" IMAGE_DIR="predictions/detection/rf_classification/647785" MASK_DIR="fov_masks" OVERWRITE=0 THRESH=0.5 MODEL_NAME="rf" MAX_SIZE=1000 qsub -V -t 1 -l twoday matlab_code/trunk/hydra/cuc/predict_apex_offsets_set.sh





