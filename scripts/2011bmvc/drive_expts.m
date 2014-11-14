% script to collate outputs from all experiments on DRIVE data
clear; clc;

decomps = {'Monogenic','2nd deriv','Haar-like','DT-CWT'};
reg_labels = {'Analytic','LinReg','Boost ','Forest'};
linestyles = {':','--','-.','-'}; % same number as #regressors
colororder = [ 0 0 1; 1 0 0; 0 0.5 0; 0 0.75 0.75 ];

outpath = 'S:\projects\mammography\matlab\papers\2011bmvc\';

% datasets = {'mammography','retinogram'};
% datasets = {'mammography'};
datasets = {'retinogram'};
for i_dataset = 1:length(datasets)
	dataset = datasets{i_dataset};

	figure(1); clf; hold on;
		plot([0,90],0.5*[1,1],'-','color',0.9*[1,1,1]);
		plot([0,90],0.9*[1,1],'-','color',0.9*[1,1,1]);
	h = []; legstrs = {};
		
	switch dataset
		case 'mammography',
			title('Synthetic Mammogram-like Images')
			resroot		= [asymmetryroot,'data\synthetic_lines\real512\results\'];
			predroot	= [asymmetryroot,'data\models\mammograms\'];
			figpath		= 'mammo';
			%		mono		g2d/clover,	haar,		dtcwt
			jobs = {'mono',		'clover',	'',			''; % analytic
					'313850',	'*313801',	'313759',	'*313751'; % linreg
					'313851',	'*313765',	'313763',	'*313764'; % boost
					'10766',	'*287361',	'313715',	'*286712'}; % forest
								% ^ or 10762
			axlim = [0,90,0,1];

		case 'retinogram',
			title('Retinogram Images')
			resroot		= [asymmetryroot,'data\retinograms\drive\test\images_extended\results\'];
			predroot	= [asymmetryroot,'data\models\vessel\'];
			figpath		= 'retina';
			%		mono		g2d/clover,	haar,		dtcwt
			jobs = {'mono',		'clover',	'',			''; % analytic
					'313836',	'*313803',	'*313804',	'*313814'; % linreg
					'313838',	'313829',	'313820',	'313813'; % boost
					'313846',	'*clover_rf','*haar_rf','*dt_rf'}; % forest
			axlim = [0,45,0,1];
	end

	fid = fopen([outpath,sprintf('%s_table.txt',dataset)],'w');

	for i_reg = 1:length(reg_labels)
		percentiles = [];
		fprintf(fid,'%s\t',reg_labels{i_reg});
		reg_jobs = jobs(i_reg,:);

		for i_decomp = 1:length(decomps)
			if ~isempty(reg_jobs{i_decomp})
				f_plot = (reg_jobs{i_decomp}(1)=='*');
				if f_plot, reg_jobs{i_decomp} = reg_jobs{i_decomp}(2:end); end
				load([resroot,reg_jobs{i_decomp},'/errors.mat']);
				if f_plot
					percentiles(end+1,:) = es.abs_percentiles;
					legstrs{end+1} = sprintf('%s %s',reg_labels{i_reg},decomps{i_decomp});
				end
				fprintf(fid,'& %6.1f\t',es.abs_median);
			else
	% 			percentiles(end+1,:) = nan(1,100);
	% 			legstrs{end+1} = sprintf('%s %s',reg_labels{i_reg},decomps{i_decomp});
				fprintf(fid,'& %6s\t','-');
			end
		end
		fprintf(fid,'\\\\\n');
		
		if ~isempty(percentiles)
			set(gca,'linestyleorder',linestyles{i_reg},...
					'colororder',colororder);
			h = [h; plot(percentiles,0.01:0.01:1)];
		end
	end % for i_reg
	fclose(fid);

	axis(axlim);
	legend(h,legstrs,'location','southeast');
	xlabel('Orientation error (degrees)'); ylabel('Cum. Freq');

	graph(1); exportfig([outpath,sprintf('figs/%s/%s_expt',figpath,dataset)],...
						{'fig','eps','pdf'});
end % for i_dataset	
	
	