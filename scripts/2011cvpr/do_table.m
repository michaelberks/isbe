% script to generate table for CVPR paper

n_sets = 1; % number of bootstrap samples (1 = no bootstrap)
n_per_set = inf; % number of samples per bootstrap (inf = all)

%% read in error data
if ~exist('error_struct','var')
dataroot = [asymmetryroot('shared'),'data/retinograms/DRIVE/test/predictions/'];

for pred_type = {'analytic','linear_regression','boosted_regression','rf'}
	for window_size = [1,3]
		for decomp = {'g1d','g2d','mono','dt'}
			if strcmp(pred_type{1},'analytic')
				datadir = [dataroot,decomp{1},'/analytic/errors/'];
			else
				datadir = [dataroot,decomp{1},'/',pred_type{1},'_',num2str(window_size),'/errors/'];
			end
			filename = [datadir,'orientation_errors.mat'];
			
			all_vec = NaN;
			centre_vec = NaN;
			if exist(filename,'file')
				load(filename);

				if (n_sets==1)
					all_vec = median(abs(prediction_errs));
					centre_vec = median(abs(centre_errs));
				else
					all_vec = zeros(1,n_sets);
					centre_vec = zeros(1,n_sets);
					
					for i = 1:n_sets
						n_samples = min(n_per_set,length(prediction_errs));
						inds = ceil(rand(1,n_samples)*length(prediction_errs));
						all_vec(i) = median(abs(prediction_errs(inds)));

						n_samples = min(n_per_set,length(centre_errs));
						inds = ceil(rand(1,n_samples)*length(centre_errs));
						centre_vec(i) = median(abs(centre_errs(inds)));
					end
				end
			end
			
			error_struct.(decomp{1}).(pred_type{1}).(['w',num2str(window_size)]) = ...
				struct('centre', centre_vec, 'all', all_vec);
		end
	end
end
end

%% dump it into a table
pred_types = {'analytic','linear_regression','boosted_regression','rf'};
decomp_types = {'g1d','g2d','mono','dt'};

pred_labels = {'Analytic','LinReg','Boost','Forest'};
% decomp_labels = {'$1^{st}$ deriv.','$2^{nd}$ deriv.','Mono.','DT-RF'};
decomp_labels = {'$G''$','$G''''$','Mono.','\dtcwt'};

filename = [asymmetryroot('shared'),'data\retinograms\DRIVE\test\predictions\retinogram_table.txt'];
fid = fopen(filename,'w');
	
fprintf(fid,'%s\n', '\begin{tabular}{l c c c c c}');
fprintf(fid,'%s\n', '\toprule');
fprintf(fid,'%s\n', '& & \multicolumn{4}{c}{Feature Type} \\');
fprintf(fid,'%s', '& ');
for i = 1:length(decomp_labels)
	fprintf(fid,'& %s ', decomp_labels{i});
end
fprintf(fid,'%s\n', '\\');
fprintf(fid,'%s\n', '%');

for i_pred_type = 1:length(pred_types)
	pred_type = pred_types{i_pred_type};
	if strcmp(pred_type,'analytic'),	n_rows = 2;
	else								n_rows = 4;
	end
	
	fprintf(fid,'%s\n', '\cmidrule{2-6}');
	s = sprintf('\\multirow{%i}{*}{%s} &',n_rows,pred_labels{i_pred_type});
	fprintf(fid,'%s\n',s);

	for window_size = [1,3]
		if strcmp(pred_type,'analytic') && window_size==3
			continue;
		end
		
		if (window_size>1)
			fprintf(fid,'  &\n');
		end
		s = sprintf('\\multirow{2}{*}{%i{$\\times$}%i} &',window_size,window_size);
		fprintf(fid,'%s\n',s);

		% print error over all points
		fprintf(fid,'    ');
		for i_decomp = 1:length(decomp_types)
			% get errors for this condition
			decomp = decomp_types{i_decomp};
			es = error_struct.(decomp).(pred_type).(['w',num2str(window_size)]);
			
			errstr = sprintf('%.1f ',mean(es.all)*180/pi);
			fprintf(fid,'%8s ',errstr);
% 			varstr = sprintf('%.1f ',var(es.all)*180/pi);
% 			fprintf(fid,'[%8s] ',varstr);
			if (i_decomp<length(decomp_types)),	fprintf(fid,'\t& ');
			else								fprintf(fid,'\\\\ \n');
			end
		end

		% print error over centre points
		fprintf(fid,'& & ');
		for i_decomp = 1:length(decomp_types)
			% get errors for this condition
			decomp = decomp_types{i_decomp};
			es = error_struct.(decomp).(pred_type).(['w',num2str(window_size)]);
			
			errstr = sprintf('(%.1f)',mean(es.centre)*180/pi);
			fprintf(fid,'%8s ',errstr);
% 			varstr = sprintf('(%.1f)',var(es.centre)*180/pi);
% 			fprintf(fid,'[%8s] ',varstr);
			if (i_decomp<length(decomp_types)),	fprintf(fid,'\t& ');
			else								fprintf(fid,'\\\\ \n');
			end
		end
		
		fprintf(fid,'%s\n', '%');
	end
end
fprintf(fid,'%s\n', '\bottomrule');
fprintf(fid,'%s\n', '\noalign{\smallskip}');
fprintf(fid,'%s\n', '\end{tabular}');

fclose(fid);

