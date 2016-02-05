function score_CompressedBreast_Example_VDM_Sample_PC()
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
rng('shuffle');
generator_status=rng;
initpath=cd(); %PC fixme
warning off all;
%store the paths for the various images necessary for the experiments
%this is the path of the Test Set
Volpara_Density_Maps_Path=strcat(initpath,'/TestSet/');
%the results of the RF models
RF_VDM_Path=strcat(initpath,'/Outputs/RF_VDM/');
%the results of the scores
SCORES_VDM_Path=strcat(initpath,'/Outputs/SCORES_VDM/');

% assign the decomposition arguments for DT_CWT
decomposition_args.decomp_type = 'dt';      %Use the dual-tree
decomposition_args.win_size = 1;            %Window size about pixel to sample features
decomposition_args.levels = 1:6;            %Levels to include
decomposition_args.feature_shape = 'rect';  %Keep this as rect, although there are other things you can play with (see sample_dt_data)
decomposition_args.feature_type = 'conj';   %Again, 'conj' works best, but can use 'mag' for magnitude only, or 'phase' (see convert_complex_representation)
decomposition_args.do_max = 0;              %Keep all 6 sub-bands (0) or only keep the maximum (1)
decomposition_args.rotate = 0;              %Try and make rotation invariant (best left switched off)
decomposition_args.use_nag = 0;             %Use the NAG toolbox to interpolate (best left swicthed off)
decomposition_args.normalise = 0;

%the number of samples that will be used for every image
% % sample_size_bank=[1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 11000 12000 13000 14000 15000 16000 17000 18000 19000 20000 21000 22000 23000 24000 25000 26000 27000 28000 29000 30000];
% sample_size_bank=[9600 13000 14000 15000 16000 17000 18000 19000 20000 21000 22000 23000 24000 25000 26000 27000 28000 29000 30000];
% sge_idx=str2num(getenv('CLUSTER_INDEX'));
% [chunk_index sample_size_indx]=ind2sub([10 30], sge_idx);
% number_of_samples=sample_size_bank(sample_size_indx);

%PC fixme
number_of_samples=1000;
PixelsScored=5000;

% load the folder list
% cd '/scratch/mopcgem2/';
% load 'D84ImageList.mat'; % it load a 14159 x 1 cell with the folder paths
cd(Volpara_Density_Maps_Path);
InitialImageList=Record_Image_Names_and_Paths('.pgm',0);
D84ImageList=unique(InitialImageList(:,2));
disp('ImageList loaded');   

%move to the path of the RF Model
cd(RF_VDM_Path);
%load the RF model (the variable that is loaded is: RF_model)
tic; load 'RF_model_VDM_nos_1000_Example'; toc;
disp('RF model loaded');           

 for i=1:length(D84ImageList)
        
        cd (D84ImageList{i,1}); %move to the appropriate subdirectory
        newpath=cd();          %get the full path
        clear ListOfImages ListOfMasks;
        ListOfImages=ReadContentOfDirectory(newpath,'densityMap'); %get a list of the VDM images
        ListOfMasks=ReadContentOfDirectory(newpath,'Segmentation'); %get a list of the VOLPARA Masks
        %loop through the images of the directory and find which ones are
        %the correct ones
        clear ind_vec;
        counter=0;
        for cntr=1:length(ListOfImages)
            if isempty(strfind(ListOfImages{cntr},'MLO'))==0
               counter=counter+1;
               ind_vec(counter)=cntr; 
            end
        end
        
        %now process the correct ones
        
        for int_cntr=1:length(ind_vec)
          % get the foreground mask
          newMask=imread(ListOfMasks{ind_vec(int_cntr)});      % read the volpara mask 
          %get the glandular and fatty tissue
          ind4=find(newMask(:)==4);
          Mask=zeros(size(newMask));
          Mask(ind4(:))=1;

          % final pixel candidates
          [ycord xcord]=find(Mask==1);

          % load the images
          img=double(imread(ListOfImages{ind_vec(int_cntr)}));     
          img=(img./65535 - 0.5) .* 2.0; %transformation to become actual density maps 
            if size(newMask)~=size(img)
               error('the mask does not match the size of the image');
            end
          %produce the final name
          name=ListOfImages{ind_vec(int_cntr)};
          position=strfind(name,'-RAW');
          name(position:end)=[];
          maskName=name;

          [ydim xdim]=size(img);

          % now compute the DT-CWT responses
          [responses] = compute_filter_responses(img, decomposition_args);

          % sample the points
          size_of_population = length(xcord);
          sample_indices = randperm(size_of_population, PixelsScored);
          sampled_features = zeros(PixelsScored,72);      %   initialise matrix of final sampled features
          sampled_features(1:length(sample_indices),:)=sample_image_features(responses, ycord(sample_indices(1:end)), xcord(sample_indices(1:end)), decomposition_args);

          %predict according to the RF_model
          test_options.predict_all = 1;
          display('Fortran prediction');
          tic;
          [y, votes, prediction_per_tree] = classRF_predict(sampled_features, RF_model, test_options);
          toc;
          
          rf_mat = convert_fortran_rf(RF_model);
          
          display('Matlab prediction')
          tic;
          [y_mat votes_mat prediction_per_tree_mat] = random_forest_class_predict(rf_mat, sampled_features, [], 0);
          toc;
          
          mismatch_class = y_mat(:) ~= y;
          mismatch_votes = votes_mat ~= votes;          
          [~, pptm] = max(prediction_per_tree_mat,[],2);
          mismatch_trees = prediction_per_tree ~= squeeze(pptm);
          
          display(['Number of mis-matched classes: ' num2str(sum(mismatch_class))]);
          display(['Number of mis-matched votes: ' num2str(sum(mismatch_votes(:,1)))]);
          display(['Number of mis-matched trees: ' num2str(sum(mismatch_trees(:)))]);
          
          RF_votes=votes;
          clear Y_hat prediction_per_tree;         

          %the votes need be normalised (200 trees) to represent
          %probability. Therefore they are divided by 2 (for 300 trees they'd be divided by 3)
          score_class_1(:,1)=sum(RF_votes(:,1)./2)/(length(RF_votes));
          score_class_2(:,1)=sum(RF_votes(:,2)./2)/(length(RF_votes));
          clear votes RF_votes;

          %move to the appropriate directory 
          cd(SCORES_VDM_Path);
          %save scores in a mat file
          save([maskName '_VDM_SCORE_' int2str(PixelsScored) '_Pixels_nos_' int2str(number_of_samples) '_Example_Model.mat'], 'score_class_1', 'score_class_2', 'generator_status');       
          clear score_class_1 score_class_2 ;   

          clear votes RF_votes;

          cd (D84ImageList{i,1}); %move back to the appropriate subdirectory
        
        end     
 end
