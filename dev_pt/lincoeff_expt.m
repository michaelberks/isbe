clc; %clear all;
figroot = [asymmetryroot,'figs/'];
figpath = figroot;

% make figure dir if not already there
if ~exist(figpath,'dir')
	mkdir(figpath);
end

bg = zeros(256,256);

sampling_args = struct(...
	'win_size',1);

% generate a synthetic line around 360deg
if ~exist('coeffs','var')
	coeffs = zeros(180,60);
	for theta = 1:180
		line_args = struct(...
			'orientation_range',theta*[1,1],...
			'width_range',8*[1,1],...
			'contrast_range',0.5*[1,1],...
			'squash_range',0*[1,1],...
			'line_type','sin',...
			'decay_rate',0);
		img = generate_line_image(bg,line_args);
		if 0
			figure(1); clf; hold on; colormap(gray(256));
				imagesc(img);
				plot(128,128,'b.');
				axis('image','ij');
		end

		dt = compute_dual_tree(img,5,0);
		coeffs(theta,:) = sample_dt_data(dt,129,128,sampling_args);
	end
end

figure(1); clf;
for c = 1:30
	subplot(10,6,c);
	plot(1:180,coeffs(:,c));
	xlim([1,180]);
end
for c = 30+(1:30)
	subplot(10,6,c);
	plot(1:180,coeffs(:,c));
	xlim([1,180]);
end
% for c = 30+(1:24)
% 	subplot(10,6,c);
% 	plot(1:180,coeffs(:,c+6)-coeffs(:,c));
% 	xlim([1,180]);
% end

return
	

%% build regressors
fprintf('Building...');
for t = 1:n_trees
	load([forest.tree_root,forest.tree_dir,sprintf('traindata%04d.mat',t)]);
	X(:,31:end) = -X(:,31:end);
	
% 	% boosted regressor
% 	fprintf('boosted...');
% 	boostreg = pt_boosted_reg_train(X,y,...
% 		'shrinkage',0.25,'n_levels',100,'n_bins',24);
	
	% linear regressor
	fprintf('linear...');
	linreg = pt_lin_reg_train(X,y);
end
fprintf('done\n');

% if exist('boostreg','var')
% 	inps = [boostreg.levels(:).input];
% 	figure(3); hist(inps,1:size(X,2));
% end

D = size(X);
x = (60*((1:60)-0.5)) * pi/180;
lines_x = ones(2,1)*(6:6:54);
figure(1); clf; 
subplot(2,1,1); hold on;
	plot((1:60)-0.5,cos(x)/16,'r.:');
	plot((1:60)-0.5,real(linreg.beta(2:end)),'r-');
	lines_y = ylim'*ones(1,length(lines_x));
	plot((1:60)-0.5,cos(x)/16,'r.:');
	plot(lines_x,lines_y,'-','color',0.8*[1,1,1]);
	plot((1:60)-0.5,real(linreg.beta(2:end)),'r-');
	xlabel('Feature'); ylabel('real(beta)'); title('Real parts of beta');
subplot(2,1,2); hold on;
	plot((1:60)-0.5,sin(x)/16,'b.:');
	plot((1:60)-0.5,imag(linreg.beta(2:end)),'b-');
	lines_y = ylim'*ones(1,length(lines_x));
	plot((1:60)-0.5,sin(x)/16,'b.:');
	plot(lines_x,lines_y,'-','color',0.8*[1,1,1]);
	plot((1:60)-0.5,imag(linreg.beta(2:end)),'b-');
	xlabel('Feature'); ylabel('imag(beta)'); title('Imaginary parts of beta');
graph(1); set(1,'paperposition',[0 0 20 15]); 
% exportfig([figpath,'linreg_coeffs']);

return

%% predict training data (get training errors)
fprintf('Predicting (train)...');
fprintf('linear...');	y_linear = linear_regressor_predict(linreg,X);
fprintf('boosted...');	y_boosted = boosted_regressor_predict(boostreg,X);
fprintf('forest...');	y_tree = mb_random_forest_reg_predict(forest,X);
fprintf('done\n');

[ignore,es_linear] = ori_error(y,y_linear); display(es_linear);
[ignore,es_boosted] = ori_error(y,y_boosted); display(es_boosted);
[ignore,es_tree] = ori_error(y,y_tree); display(es_tree);


%% generate test data
if (1)
	if ~exist('X_test','var')
		% sample test data myself if it doesn't exist already
		fprintf('Generating test data...');
		sampling_args = rmfield(sampling_args,'detection_type');
		sampling_args.num_samples = 4000;
		sampling_args.task_id = 1;
		[X_test,y_test] = sample_saved_dt_line_data(sampling_args);
		fprintf('done\n');
	end
else
	% use Mike's sampled data
	X_test = [];
end

%% predict test data (get test error)
fprintf('Predicting (test)...');
y_linear2 = []; y_linear = []; y_boosted = []; y_tree = [];
if isempty(X_test)
	% use Mike's sampled data
	
	y_test = []; % all test targets in one place
	datapath = 'A:\data\synthetic_data\real512_dt\test\';
	datadir = dir([datapath,'X_*.mat']);
	tb = timebar('title','Processing datasets','limit',length(datadir));
	for i = 1:length(datadir)
		% load test data (X and y)
		X_test = u_load([datapath,datadir(i).name]);
		X_test = convert_dt_representation(X_test(:,5,:,:),...
					'feature_type',sampling_args.feature_type,...
					'win_size',sampling_args.win_size);
				
		% this data has outputs in range [0..pi]
		% double to [-pi..pi]
		y_tmp = u_load([datapath,strrep(datadir(i).name,'X_','y_')]);
		y_test = [y_test; exp(sqrt(-1)*angle(y_tmp)*2)];

% 		% predict with tree
% 		y_tree = [y_tree; mb_random_forest_reg_predict(forest,X_test)];

		% predict with boosted regressor
		y_boosted = [y_boosted; boosted_regressor_predict(boostreg,X_test)];
		
		timebar(tb,'advance');
	end
	timebar(tb,'close'); clear tb;
	
	% set error file suffix
	err_suffix = '_mb';
else
	% use my sampled data
	
	fprintf('linear...');	y_linear = [y_linear; linear_regressor_predict(linreg,X_test)];
	fprintf('boosted...');	y_boosted = [y_boosted; boosted_regressor_predict(boostreg,X_test)];
	fprintf('forest...');	y_tree = [y_tree; mb_random_forest_reg_predict(forest,X_test)];

	% set error file suffix
	err_suffix = '';
end
fprintf('done\n');

% save the errors themselves
save([figpath,'outputs',err_suffix,'.mat'],...
	 'y_test','y_linear','y_boosted','y_tree');

[e_linear,err_stats_linear] = ori_error(y_test,y_linear);
[e_boosted,err_stats_boosted] = ori_error(y_test,y_boosted);
[e_tree,err_stats_tree] = ori_error(y_test,y_tree);

%% output

% display on screen
display(err_stats_linear);
display(err_stats_boosted);
display(err_stats_tree);

% dump errors to log file
filename = [figpath,'errors',err_suffix,'.txt'];
fid = fopen(filename,'w');
	fprintf(fid,'User = %s\n',get_username());
	fprintf(fid,'%s\n',evalc('display(err_stats_linear);'));
	fprintf(fid,'%s\n',evalc('display(err_stats_boosted);'));
	fprintf(fid,'%s\n',evalc('display(err_stats_tree);'));
fclose(fid);

figure(1); clf;
	plot(err_stats_linear.abs_percentiles,(1:100)/100,'b-',...
		 err_stats_boosted.abs_percentiles,(1:100)/100,'r-',...
		 err_stats_tree.abs_percentiles,(1:100)/100,'g-');
	axis([0,90,0,1]);
legend({'linear','boosted','forest (40x16k)'},'location','southeast');
graph(1); exportfig([figpath,'cum_freq',err_suffix]);

return

% plot betas (good and bad)
figure(2); clf;
	plot(abs(beta));
graph(2); exportfig([figpath,'betas_good']);
	plot(abs(beta_bad));
graph(2); exportfig([figpath,'betas_bad']);

% get histogram over split variables and plot
n_trees = 1;
hists = zeros(size(X_test,2),n_trees);
for t = 1:n_trees
    tree = u_load([forest.tree_root,forest.tree_dir,forest.trees{t}]);
	hists(:,t) = hist(tree.var(tree.var>0),1:60);
	hists(:,t) = hists(:,t)/sum(hists(:,t));
end
figure(2); clf;
	plot(hists);
graph(2); exportfig([figpath,'betas_treehists']);

profstat = profile('status');
if strcmp(profstat.ProfilerStatus,'on')
	profile off; profile report;
end
