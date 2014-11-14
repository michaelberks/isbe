function seed_randomizers(args)

% set random seed based on task id and clock to ensure unique trees created
if isfield(args,'rand_seed') && ~isempty(args.rand_seed);
    rand('twister', args.rand_seed);
	randn('state', args.rand_seed);
else
    if ~isfield(args,'task_id')
        error('args must have either a task_id or non-empty rand_seed field');
    end
    
    rand('twister', sum(args.task_id*clock)); %#ok
	randn('state', sum(args.task_id*clock));
end
