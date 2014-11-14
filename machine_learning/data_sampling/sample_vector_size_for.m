function sz = sample_vector_size_for(args,sampling_args)

sz = NaN;

check_decomp_type(args.decomp_type);

switch args.decomp_type
	case 'dt'
		sz = args.num_levels*compute_dt_feature_size(sampling_args);
	case 'mono'
		sz = 3*args.num_levels*args.win_size^2;
	case {'g2d','clover'}
		sz = 3*length(args.sigma_range)*args.win_size^2;
	case {'g1d','haar'}
		sz = 2*length(args.sigma_range)*args.win_size^2;
	case 'linop'
		sz = (args.num_angles - args.do_max*(args.num_angles-1))*...
			 (args.num_levels*args.win_size^2);
	case 'pixel'
		sz = args.win_size^2;
	otherwise
		error(['Sample vector size not defined for ',args.decomp_type]);
end    
