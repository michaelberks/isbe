function [model_definition] = get_gm_tex_model_definition(in_string)
%
% GET_GM_TEX_MODEL_DEFINITION
%
% This is the definition of a model building struct. It is used to verify that
% the model building struct used to define the model meets the correct
% specification.
%
% No paramters can be set in this function; parameters that cannot be set by
% the user (such as the model version) are set in build_gm_texture_model.
%
% There are two sections to the model definition: Parameters that the user must set
% and parameters that the user must not set (i.e. they are set in code).
%
% get_gm_tex_model_definition('disp');
%
% displays the parameters the user must define (use the ; to suppress the whole
% list of parameters)
%
% get_gm_tex_model_definition('disp_all');
%
% DEVELOPERS ONLY: displays the all parameters that must be defined (use
%	the ; to suppress the whole list of parameters)

user_params = cell(0); % the params that the user must set
non_user_params = cell(0); % the params that the user must not set
user_params_docs = cell(0); % some documentation for each user parameter
non_user_params_docs = cell(0); % some documentation for each non-user parameter

% Required user model parameters
% -------------------------------------------

% Clustering Parameters
user_params{end+1} = 'k';
user_params_docs{end+1} 				= 'The number of model components in the model';
user_params{end+1} = 'NextDataFunction';
user_params_docs{end+1} 				= 'The name of the function to get the next data for clustering';
user_params{end+1} = 'NextDataFunctionArgs.ImageDir';
user_params_docs{end+1} 				= 'The directory containing the training images';
user_params{end+1} = 'NextDataFunctionArgs.TempStorageDir';
user_params_docs{end+1} 				= 'An empty directory that can be used for temporary storage (you can reuse temp dirs)';
user_params{end+1} = 'NextDataFunctionArgs.WindowSize';
user_params_docs{end+1} 				= 'The size of the window to be modelled -- must be odd';
user_params{end+1} = 'NextDataFunctionArgs.MaxMemory';
user_params_docs{end+1} 				= 'The maximum amount of memory to use for the next chunk of data';
user_params{end+1} = 'ClusteringFunction';
user_params_docs{end+1} 				= 'The name of the clustering function';
user_params{end+1} = 'TempDirPath';
user_params_docs{end+1} 				= 'A directory that can be used for temporary storage (you can reuse temp dirs)';
user_params{end+1} = 'MaxMemoryFinalClusteringInMB';
user_params_docs{end+1} 				= 'The maximum amount of RAM to be used in the final clustering';
user_params{end+1} = 'ClusteringFunctionArgs.Iterative';
user_params_docs{end+1} 				= 'Determines if the iterative or non-iterative variant of the clustering function is used';
user_params{end+1} = 'ClusteringFunctionArgs.SmallClusterSize';
user_params_docs{end+1} 				= 'The size of small clusters (small clusters are removed)';
user_params{end+1} = 'ClusteringFunctionArgs.PropRepresentativePoints';
user_params_docs{end+1} 				= 'The proportion of each cluster to return as representative points';
user_params{end+1} = 'FinalPassClusteringFunction';
user_params_docs{end+1} 				= 'The name of the clustering function to be used in the final clustering';
user_params{end+1} = 'FinalPassClustFunctionArgs.Iterative';
user_params_docs{end+1} 				= 'Determines if the iterative or non-iterative variant of the clustering function is used';
user_params{end+1} = 'FinalPassClustFunctionArgs.PropRepresentativePoints';
user_params_docs{end+1} 				= 'The proportion of each cluster to return as representative points';
user_params{end+1} = 'FinalPassClustFunctionArgs.SmallClusterSize';
user_params_docs{end+1} 				= 'The size of small clusters (small clusters are removed)';
user_params{end+1} = 'MaxRunTime';
user_params_docs{end+1} 				= 'The maximum amount of time to the clustering for before moving on to the final clustering';

% Model PCA parameters
user_params{end+1} = 'PCA.PerformPCA';
user_params_docs{end+1} 				= 'Whether or not to perform a PCA (set to true or false); if false, other PCA parameters must be set to empty, and will be ignored';
user_params{end+1} = 'PCA.Keep';
user_params_docs{end+1} 				= 'The proportion of the variance of the total variance to keep (must be a value in the interval [0,1])';
user_params{end+1} = 'PCA.MaxRunTime';
user_params_docs{end+1} 				= 'The maximum amount of time (in hours) to spend estimating the PCA';

% Model Description Parameters
user_params{end+1} = 'Description';
user_params_docs{end+1} 				= 'A short description of the model';

% Model Save Parameters
user_params{end+1} = 'ModelSavePath';
user_params_docs{end+1} 				= 'The full path to the location to save the model in';

% Required non-user parameters
% -------------------------------------------
non_user_params{end+1} = 'ModelVersion';
non_user_params_docs{end+1} 				= 'The version number of the model';
non_user_params{end+1} = 'ModelBuiltOn';
non_user_params_docs{end+1} 				= 'The date the model was built on';

% Define the model:
model_definition = {user_params{:} non_user_params{:}};

if nargin == 1
	if strcmp(in_string, 'disp') | strcmp(in_string, 'disp_all')
        disp(' ')
		disp('The parameters the user must define in the model building struct are:')
        disp('---------------------------------------------------------------------')
		for i = 1 : length(user_params)
			disp(user_params{i})
			disp(['    - ' user_params_docs{i}])
			disp(' ')
		end
	end
	if strcmp(in_string, 'disp_all')
		disp('The parameters the developers must define in the code are:')
        disp('----------------------------------------------------------')
		for i = 1 : length(non_user_params)
			disp(non_user_params{i})
			disp(['    - ' non_user_params_docs{i}])
			disp(' ')
		end
	end
end
