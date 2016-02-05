function Run_Compute_Train_Score_Example_Experiment_VDM_PC(data_path)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%data_path=cd();
tic; compute_CompressedBreast_DT_CWT_Example_VDM_PC(data_path); toc;
disp('finished computing the DT_CWT coefficients of the training set images');
cd(data_path);
tic; train_CompressedBreast_RF_model_Example_VDM_PC(); toc;
disp('finished training the RF model');
cd(data_path);
tic; score_CompressedBreast_Example_VDM_Sample_PC; toc;
disp('finished scoring the test images');
cd(data_path);