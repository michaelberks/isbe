function model = mb_build_gaussian_pyramid_model(varargin)

args = u_packargs(varargin,... % the user's input
			 '0', ... % strict mode
			 {'LoadIdx',...
             'Idx',...
             'ImageDir',...
             'CutOffLevel'},...
             'ImageList', [],...
             'WindowSize', 15,...
             'NumLevels', 5,...
             'NumOrientations', 5,...
             'MaxMemory', 128,...
             'SaveFile', []);

clear varargin;

model.Means = cell(args.CutOffLevel, 1);
model.Covars = cell(args.CutOffLevel, 1);

%for each level
for level = 1:args.CutOffLevel
    
    %Get data
    data_args.LoadIdx = args.LoadIdx;
    data_args.Idx = args.Idx;
    data_args.ImageDir = args.ImageDir;
    data_args.Level = level;
    data_args.MaxMemory = args.MaxMemory;
    data = get_complete_data_for_clustering(data_args);
    
    %quick hack
    data(sum(isnan(data), 2) > 0, :) = []; 
    sum(isnan(data(:)))
    
    % Standardise data
    [data, means, stds] = st_standardise_data(data);

    %Get mean and covariance of data
    model.Means{level} = mean(data);
    model.Covars{level} = cov(data);
    model.St_Means{level} = means;
    model.St_Stds{level} = stds;
    
end

if ~isempty(args.SaveFile)
    save(args.SaveFile, 'model');
end
    