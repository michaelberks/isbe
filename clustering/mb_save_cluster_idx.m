function [ClusterIdxArgs] = mb_save_cluster_idx(varargin)


args = u_packargs(varargin, ... % what the user sent us
			'0',...
			{...  % The mandatory arguments
            'NumberSets',...
            'NextIdxFunctionArgs',...
            'NextIdxFunction'},... % The optional arguments
            'MinimumPoints', 2000,...
            'IdxFile', 'cluster_idx',...
            'IdxDir', [mberksroot, 'background/idx/']);
clear varargin;

if ~isdir(args.IdxDir)
    mkdir(args.IdxDir)
end
if args.IdxDir(end) ~= filesep
	args.IdxDir = [args.IdxDir filesep];
end

number_sets = args.NumberSets;

% pre-allocate space for set filenames
ClusterIdxArgs.Files(number_sets).name = [];

retained_sets = 0;
for ii = 1:number_sets
    
    [cluster_idx, finished_sampling_flag] = feval(args.NextIdxFunction, args.NextIdxFunctionArgs);
    if sum(~isnan(cluster_idx(:))) > args.MinimumPoints
        % update the count of retained sets
        retained_sets = retained_sets + 1;
        
        %if we've got enough points save the indices
        idx_name = [args.IdxDir, args.IdxFile, zerostr(retained_sets, 3)];
        save(idx_name, 'cluster_idx');
        
        %record the set name
        ClusterIdxArgs.Files(retained_sets).name = idx_name;
        
    else
        % If this is too small, so are all future sets, so break now
        break;
    end
    if finished_sampling_flag 
        %if there is no more data to sample, then quit this loop
		break;
    end
end
%make next idx function clear up any temporary files
args.NextIdxFunctionArgs.CleanUp = (1==1);
feval(args.NextIdxFunction, args.NextIdxFunctionArgs); % perform the clean up

%return the number of sets saved and remove any blanks from the end of
%files names
ClusterIdxArgs.NumberSets = retained_sets;
ClusterIdxArgs.Files(retained_sets+1:end) = [];

    
    
    