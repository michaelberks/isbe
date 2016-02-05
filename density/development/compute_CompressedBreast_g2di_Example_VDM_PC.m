function compute_CompressedBreast_g2d_Example_VDM_PC(initpath)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 
%  function compute_CompressedBreast_DT_CWT_Example_VDM_PC()
% 
%  syntax: compute_CompressedBreast_DT_CWT_Example_VDM_PC()
%
% Emmanouil Moschidis 10/09/2015
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

rng('shuffle');
generator_status=rng;

% warning('off', 'ASYM:unexpectedArgument'); %switch off the warnings
warning off all;

%store the paths for the various images necessary for the experiments
Volpara_Density_Maps_Path=strcat(initpath,'/TrainSet/');
%The computed DT_CWT feature images
DT_CWT_VDM_Path=strcat(initpath,'/Outputs/DT_CWT_VDM/');
%the results of the RF models
RF_VDM_Path=strcat(initpath,'/Outputs/RF_VDM/');

% first assign the decomposition arguments
decomposition_args.decomp_type = {'g2dia', 'h2dai'}; %Downsampling/Interpolating form of Gaussian derivatives
decomposition_args.win_size = 1;            %Window size about pixel to sample features
decomposition_args.sigma_range = [1 6];     %Finest scale Gaussian sigma and number of levels
decomposition_args.num_angles = 6;
decomposition_args.do_max = 0;              %Keep all 6 sub-bands (0) or only keep the maximum (1)
decomposition_args.rotate = 0;              %Try and make rotation invariant (best left switched off)
decomposition_args.normalise = 0;
%cluster stuff
% sample_size_bank=[1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 11000 12000 13000 14000 15000 16000 17000 18000 19000 20000 21000 22000 23000 24000 25000 26000 27000 28000 29000 30000];
% sge_idx=str2num(getenv('CLUSTER_INDEX'));
% [i_indx sample_size_indx]=ind2sub([500 30], sge_idx);
% number_of_samples=sample_size_bank(sample_size_indx);
init=1;

%PC stuff
number_of_samples=1000;

 
     %we will run the experiment on the VDM images
     % get the VDM images
     cd(strcat(Volpara_Density_Maps_Path,'CANCERS/'));
     CancerDensityMapsMat=Record_Image_Names_and_Paths('densityMap',0); %get the density maps and the paths
     cd(strcat(Volpara_Density_Maps_Path,'CONTROLS/')); 
     NormalDensityMapsMat=Record_Image_Names_and_Paths('densityMap',0); %get the density maps and the paths

     %get the VOLPARA segmentation masks
     cd(strcat(Volpara_Density_Maps_Path,'CANCERS/'));
     CancerMasksMat=Record_Image_Names_and_Paths('Segmentation', 0); % get the segmentation masks names and paths
     cd(strcat(Volpara_Density_Maps_Path,'CONTROLS/'));
     NormalMasksMat=Record_Image_Names_and_Paths('Segmentation', 0); % get the segmentation masks names and paths
       
     % all the necessary images (first cancers then controls)    
     MasksMat=cat(1,CancerMasksMat,NormalMasksMat);
     DensityMapsMat=cat(1, CancerDensityMapsMat, NormalDensityMapsMat);
 
%loop through the images
     
for i=1:length(DensityMapsMat)%PC loop
%         i=i_indx; % cluster stuff
        tic;
       
        % get the foreground mask
        cd(MasksMat{i,2});                 % go to the appropriate path
        newMask=imread(MasksMat{i,1});     % load the mask
        
        % load the images
        cd(DensityMapsMat{i,2});
        img=double(imread(DensityMapsMat{i,1}));
	    % Volpara normalisation transform
        img=(double(img)./ 65535 - 0.5) .* 2.0;

%         img=imresize(img,size(newMask));
        if size(newMask)~=size(img)
           error('the mask does not match the size of the image');
        end


        %get the glandular and fatty tissue
        ind4=find(newMask(:)==4);   % fully compressed tissue
        Mask=zeros(size(newMask));
        Mask(ind4(:))=1;
   
        %final pixel candidates
        [ycord xcord]=find(Mask==1);

        %produce the final name
        name=DensityMapsMat{i,1};
        position=strfind(name,'_hint');
        name(position:end)=[];
        maskName=name;
        
%         if ~isempty(strfind(name, 'RMLO'))
%            img=fliplr(img);  
%         end
 	 
        % now compute the DT-CWT responses
        [responses] = compute_filter_responses(img, decomposition_args);

        %sample the points
        size_of_population = length(xcord);
        % sample_indices = randperm(size_of_population, number_of_samples);
	    sample_indices = randi(size_of_population,1, number_of_samples); %  allow replacement
        sampled_features = zeros(number_of_samples,72);                  %  initialise matrix of final sampled features
        sampled_features(1:length(sample_indices),:) = sample_image_features(responses, ycord(sample_indices(1:end)), xcord(sample_indices(1:end)), decomposition_args);
             
        % PC
        cd (DT_CWT_VDM_Path); 
        save([maskName '-VDM-g2di_Example_nos_' int2str(number_of_samples) '_init_' int2str(init) '.mat'], 'sampled_features','generator_status');
        clear sampled_features;
        clear sample_indices;
        clear responses img newMask;
        
        toc;
       cd(initpath); %PC fixme
end