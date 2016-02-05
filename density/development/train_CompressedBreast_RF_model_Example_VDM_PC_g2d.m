function train_CompressedBreast_RF_model_Example_VDM_PC_g2d()
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

initpath=cd(); %PC fixme

%store the paths for the various images necessary for the experiments
Volpara_Density_Maps_Path=strcat(initpath,'/TrainSet/');

%The computed DT_CWT feature images
DT_CWT_VDM_Path=strcat(initpath,'/Outputs/DT_CWT_VDM/');

%the results of the RF models
RF_VDM_Path=strcat(initpath,'/Outputs/RF_VDM/');

%this is only valid for the 1-1 dataset
cancer_images=1:1;
control_images=1:1;

%the number of samples that will be used for every image
% sample_size_bank=[1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 11000 12000 13000 14000 15000 16000 17000 18000 19000 20000 21000 22000 23000 24000 25000 26000 27000 28000 29000 30000];
% sge_idx=str2num(getenv('CLUSTER_INDEX'));
% sample_size_indx=sge_idx;
% number_of_samples=sample_size_bank(sample_size_indx);

%PC stuff
number_of_samples=1000;

init=1;

% We sample to create balanced training sets. 
train_size=1;
training_labels=ones(2*train_size*number_of_samples,1);
training_labels((train_size*number_of_samples+1):end)=2;
            
        
        %now get the cancer VDM images
        ListOfCancerImages=ReadContentOfDirectory(DT_CWT_VDM_Path, 'Cancer');       % get a list of the Cancer VDM images  
        %first get the normal VDM images
        ListOfControlImages=ReadContentOfDirectory(DT_CWT_VDM_Path, 'Control');       % get a list of the Normal VDM images

        %filter out images from potentially other experiments
        intcounter=0; 
            for i=1:length(ListOfCancerImages)
                if ~isempty(strfind(ListOfCancerImages{i},['-VDM-g2d_Example_nos_' int2str(number_of_samples) '_init_' int2str(init) '.mat'])) 
                    intcounter=intcounter+1; 
%                     new_cancer_fold_indx(intcounter)=i; 
                    new_cancer_list{intcounter}=ListOfCancerImages{i};
                end 
            end
            
            
            intcounter=0; 
            for i=1:length(ListOfControlImages)
                if ~isempty(strfind(ListOfControlImages{i},['-VDM-g2d_Example_nos_' int2str(number_of_samples) '_init_' int2str(init) '.mat'])) 
                    intcounter=intcounter+1; 
%                     new_control_ind(intcounter)=i; 
                    new_control_list{intcounter}=ListOfControlImages{i};
                end 
            end
                        
            
            % I can preallocate the main matrices used for training the RF classifier
            training_features_A=zeros(train_size*number_of_samples,72);
            training_features_B=zeros(train_size*number_of_samples,72);

            
            % we find the training set 
            cd(DT_CWT_VDM_Path);
            cancer_training_set=cancer_images;    % indices of the cancer images used for training
            control_training_set=control_images;  % indices of the control images used for training
            
            
            for counter=1:length(cancer_training_set)
                load(new_cancer_list{cancer_training_set(counter)});
                training_features_A( ((counter-1)*number_of_samples+1):(counter*number_of_samples), :)=sampled_features;
                clear sampled_features;
            end
            
            %sampling from the controls is a bit more complicated 
%             for counter=1:length(control_training_set)
%                 if counter<length(control_training_set)
%                     load(new_control_list{control_training_set(counter)});
%                     training_features_B( ((counter-1)*ceil(number_of_samples/3)+1):(counter*ceil(number_of_samples/3)), :)=sampled_features;
%                     clear sampled_features;
%                 else
%                     extra_samples=(counter* ceil(number_of_samples/3))-(train_size*number_of_samples);
%                     load(new_control_list{control_training_set(counter)});
%                     training_features_B( ((counter-1)*ceil(number_of_samples/3)+1):(train_size*number_of_samples), :)=sampled_features(1:(ceil(number_of_samples/3)-extra_samples),:);                   
%                     clear sampled_features;
%                 end
%             end
            %not realistic at this point but just for the sake of of having
            %the approprate sample from controls as well
            for counter=1:length(control_training_set)
                load(new_control_list{control_training_set(counter)});
                training_features_B( ((counter-1)*number_of_samples+1):(counter*number_of_samples), :)=sampled_features;
                clear sampled_features;
            end            

                training_features=cat(1, training_features_A, training_features_B);
                clear training_features_A training_features_B;
                
               
                %train the RF model with 200 trees (not sure what will happen in terms of memory)
                disp('Ready to train the RF model...'); 
                tic; RF_model = classRF_train(training_features, training_labels, 200); toc;
                cd(RF_VDM_Path);
                save(['RF_model_VDM_nos_' int2str(number_of_samples) '_Example.mat'], 'RF_model', '-v7.3'); 
                clear RF_model;
                

                