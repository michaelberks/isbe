function [model] = mb_build_gm_texture_model(model_building_struct)
%
% BUILD_GM_TEXTURE_MODEL Builds a Gaussian Mixture Model of Texture
%
% This function builds a Gaussian Mixture Model of texture.
%
% + model_building_struct
%		a struct containing the required information for building a texture model;
%		this is detailed in GET_GM_TEX_MODEL_DEFINITION.
%		By convention, an M-file is written to build the model_building_struct
% 		(i.e. the M-file acts as a parameters file).
%
% The function returns the model.

% Set the non-user model parameters
%-----------------------------------------------
%-----------------------------------------------

% set the model version
model_building_struct.ModelVersion = 1;

% get today's date
model_building_struct.ModelBuiltOn = date;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MB: I can't be bothered with all this checking crap, just build the model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Make sure that the model as defined by the user's model building struct meets the
% current model's specification
% -------------------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------------------
%
% get the model definition
% req_params = get_gm_tex_model_definition;
% 
% % take a copy of the struct names, we'll use this to record any params that don't exist
% non_existant_params = req_params;
% 
% % determine which of the required parameters are not defined
% cells_to_remove = [];
% for i = 1 : length(req_params)
% 	% try to evaluate the sull struct name, if an error is thrown (i.e. the full struct name can't be evaluated)
% 	% then the full struct field
% % 	struct_field_exists = (1==1); % assume that it does exist
% % 	eval(['model_building_struct.' req_params{i} ';'], 'struct_field_exists = (1==0);');
% % 	if struct_field_exists
% % 		% mark this struct name to be recorded as existing
% % 		cells_to_remove(end+1) = i;
% %     end
%     struct_field_exists = (1==1); % assume that it does exist
% 	eval(['model_building_struct.' req_params{i} ';'], 'struct_field_exists = (1==0);');
%     if struct_field_exists
% 		% mark this struct name to be recorded as existing
% 		cells_to_remove(end+1) = i;
%     end
% end
% 
% % % remove all the existing parameters to create the non-existant list
% non_existant_params(cells_to_remove) = [];
% 
% % see if we have any non-existant parameters
% if ~isempty(non_existant_params)
% 	disp('The following parameters were not defined in the model definition:')
% 	for i = 1: length(non_existant_params)
% 		disp(non_existant_params{i});
% 	end
% 	error('There were undefined parameters');
% end
% 
% % Build the model
% %--------------------------------------------------
% %--------------------------------------------------
% 
% % determine if we need to estimate the global PCA
% pca_info_struct = [];
% if model_building_struct.PCA.PerformPCA
% 	% estimate the global PCA
% 	[pca_info_struct] = gm_estimate_global_pca('ImageDir', model_building_struct.NextDataFunctionArgs.ImageDir,...
% 															'TempStorageDir', model_building_struct.NextDataFunctionArgs.TempStorageDir,...
% 															'WindowSize', model_building_struct.NextDataFunctionArgs.WindowSize,...
% 															'MaxRunTime', model_building_struct.PCA.MaxRunTime);
% 	
% 	% now work out the eigenvectors and eigenvalues for the value of Keep selected
% 	pca_info_struct.Keep = model_building_struct.PCA.Keep;
%     pca_info_struct.NModes = model_building_struct.PCA.NModes;
% 	pca_info_struct.ComputePCA = 1;
% 	pca_info_struct = st_incremental_pca(pca_info_struct);
% 	
% 	% pca_info_struct is included in model below, so that we have access to the eigenvalues and mean vector when performing synthesis
% 	
% 	% set the function and args for the function to project training data into the PCA space
% 	model_building_struct.DataTransformationFunction = 'cjr_to_pca_space';
% 	model_building_struct.DataTransformationFunctionArgs = pca_info_struct;
% 	
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% perform the clustering
model = mb_gmm_cluster_large_data_set(model_building_struct);

% % Store any PCA info in the model
% if ~isempty(pca_info_struct)
% 	model.pca_info_struct = pca_info_struct;
% end

% Store the model building information in the model
model.BuildInfo = model_building_struct;

% set the WindowSize field of the model
model.WindowSize = model_building_struct.NextDataFunctionArgs.WindowSize;

% Save the model
%--------------------------------------------------
%--------------------------------------------------

% Save the model
save(model_building_struct.ModelSavePath, 'model');

% Report back to the user
%--------------------------------------------------
%--------------------------------------------------
disp('The model with the description:')
disp(model_building_struct.Description)
disp('Has been saved to:')
disp(model_building_struct.ModelSavePath)