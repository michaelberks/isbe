function Run_Compute_Train_Score_Example_Experiment_VDM_PC_g2d(data_path)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%data_path=cd();
tic; compute_CompressedBreast_g2d_Example_VDM_PC(data_path); toc;
disp('finished computing the Gaussian derivatives coefficients of the training set images');
cd(data_path);
tic; train_CompressedBreast_RF_model_Example_VDM_PC_g2d(); toc;
disp('finished training the RF model');
cd(data_path);
tic; score_CompressedBreast_Example_VDM_Sample_PC_g2d(); toc;
disp('finished scoring the test images');
cd(data_path);