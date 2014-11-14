% Script for predicting orientation in retinal images and producing stuff
% for the BMVC paper.
%
% Note that for these experiments I'm using regressors trained with 20k
% data points, usually on a 1x1 grid though some will use 3x3

clc; %clear all;
figroot = [asymmetryroot,'figs/'];

clear tree;

% profile clear; profile on;

% load most recently created forest if none specified
forest_root		= [asymmetryroot,'data/line_orientation_rfs/'];
forest_dir		= dir([forest_root,'pc*']);
forest_job		= [forest_dir(end).name,'/'];
% forest_job		= ['pc20110307T105747/']; % 20x4000
forest_job		= ['pc20110307T145627/']; % 40x16000
% forest_job		= ['pc20110307T163301/']; % 200x8000, dabs
% forest_job		= ['pc20110308T123013/']; % 200x8000, mabs
display(forest_job);
forest_dir		= dir([forest_root,forest_job,'/random_forest*.mat']);
forest_fname	= [forest_root,forest_job,forest_dir(end).name];
forest			= u_load(forest_fname);

n_trees = length(forest.trees);
n_trees = 1;
display(n_trees);

% forest.trees = forest.trees(1:n_trees);

figpath = [forest_root,forest_job,'figs/'];
% make figure dir if not already there
if ~exist(figpath,'dir')
	mkdir(figpath);
end

% get sampling arguments
sampling_args	= u_load([forest_root,forest_job,'/sampling_args.mat']);

%% build regressors
fprintf('Building...');
for t = 1:n_trees
	load([forest.tree_root,forest.tree_dir,sprintf('traindata%04d.mat',t)]);
% 	X(:,31:end) = -X(:,31:end);
% 	X = X(:,[12+(1:6),31:end]);
% 	X = X(:,[12+(1:6),12+(31:36)]);
% 	X = X(:,[12+(1:6)]);
% 	X = X(:,[12+(31:36)]);
% 	X = X(:,[1:31]);
% 	X = X(:,[31:end]);
	
% 	% boosted regressor
% 	fprintf('boosted...');
% 	boostreg = pt_boosted_reg_train(X,y,...
% 		'shrinkage',0.25,'n_levels',100,'n_bins',24);
	
	% linear regressor
	fprintf('linear...');
	linreg = pt_lin_reg_train(X,y);
	
% 	% logistic regressor
% 	fprintf('logistic...');
% 	logreg_cos = pt_log_reg_train(X,cos(angle(y))*0.5+0.5);
% 	logreg_sin = pt_log_reg_train(X,sin(angle(y))*0.5+0.5);
end
fprintf('done\n');

%% predict training data (and get training errors)
fprintf('Predicting (train)...');
% fprintf('forest...');	y_tree = mb_random_forest_reg_predict(forest,X);
fprintf('linear...');	y_linear = linear_regressor_predict(linreg,X);
% fprintf('logistic...');	y_logistic = complex(logistic_regressor_predict(logreg_cos,X)*2-1,logistic_regressor_predict(logreg_sin,X)*2-1 );
% fprintf('boosted...');	y_boosted = boosted_regressor_predict(boostreg,X);
fprintf('done\n');

[angerr,es_linear] = ori_error(y,y_linear); display(es_linear);
% [angerr,es_logistic] = ori_error(y,y_logistic); display(es_logistic);
% [ignore,es_boosted] = ori_error(y,y_boosted); display(es_boosted);
% [ignore,es_tree] = ori_error(y,y_tree); display(es_tree);

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
y_linear = []; y_logistic = []; y_boosted = []; y_tree = [];
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
% 	fprintf('forest...');	y_tree = [y_tree; mb_random_forest_reg_predict(forest,X_test)];
	fprintf('linear...');	y_linear = [y_linear; linear_regressor_predict(linreg,X_test)];
% 	fprintf('logistic...');	y_logistic = [y_logistic; complex(logistic_regressor_predict(logreg_cos,X_test)*2-1,logistic_regressor_predict(logreg_sin,X_test)*2-1 )];
% 	fprintf('boosted...');	y_boosted = [y_boosted; boosted_regressor_predict(boostreg,X_test)];

	% set error file suffix
	err_suffix = '';
end
fprintf('done\n');

% save the errors themselves
save([figpath,'outputs',err_suffix,'.mat'],...
	 'y_test','y_linear','y_logistic','y_boosted','y_tree');

%% display on screen
[e_linear,es_linear] = ori_error(y_test,y_linear); display(es_linear);
% [e_linear,es_logistic] = ori_error(y_test,y_logistic); display(es_logistic);
% [e_boosted,es_boosted] = ori_error(y_test,y_boosted); display(es_boosted);
% [e_tree,es_tree] = ori_error(y_test,y_tree); display(es_tree);

%% dump errors to log file
filename = [figpath,'errors',err_suffix,'.txt'];
fid = fopen(filename,'w');
	fprintf(fid,'User = %s\n',get_username());
	fprintf(fid,'%s\n',evalc('display(es_linear);'));
% 	fprintf(fid,'%s\n',evalc('display(es_logistic);'));
% 	fprintf(fid,'%s\n',evalc('display(es_boosted);'));
% 	fprintf(fid,'%s\n',evalc('display(es_tree);'));
fclose(fid);

%% cumulative error plots
figure(9); clf; hold on;
	plot(es_linear.abs_percentiles,(1:100)/100,'b-');
% 	plot(es_logistic.abs_percentiles,(1:100)/100,'m-');
% 	plot(es_boosted.abs_percentiles,(1:100)/100,'r-');
% 	plot(es_tree.abs_percentiles,(1:100)/100,'g-');
	axis([0,90,0,1]);
legend({'linear','logistic','boosted','forest (40x16k)'},'location','southeast');
% graph(1); exportfig([figpath,'cum_freq',err_suffix]);

%% histogram of variables selected by boosted regressor 
if exist('boostreg','var')
	inps = [boostreg.levels(:).input];
	figure(3); hist(inps,1:size(X,2));
end

%% linear regression coefficients
D = size(X,2)/2;
x = (60*((1:D)-0.5)) * pi/180;
cx = cos(x); cx = cx/(cx*cx'); sx = sin(x); sx = sx/(sx*sx');
lines_x = ones(2,1)*(6:6:D-6);
yL = 0.1*[-1,1];
figure(1); clf; 
subplot(2,2,1); hold on;
	plot((1:D)-0.5,real(linreg.beta(2:D+1)),'r-');
	plot(lines_x,yL'*ones(1,length(lines_x)),'-','color',0.8*[1,1,1]);
	plot((1:D)-0.5,cx,'r.:');
	xlabel('Scale'); ylabel('Re(\beta) = cos(2\theta)'); title('Magnitude');
	set(gca,'xtick',[3:6:D],'xticklabel',num2str((1:5)')); ylim(yL);
subplot(2,2,2); hold on;
	plot((1:D)-0.5,imag(linreg.beta(2:D+1)),'b-');
	plot(lines_x,yL'*ones(1,length(lines_x)),'-','color',0.8*[1,1,1]);
	plot((1:D)-0.5,sx,'b.:');
	xlabel('Scale'); ylabel('Im(\beta) = sin(2\theta)'); title('Magnitude');
	set(gca,'xtick',[3:6:D],'xticklabel',num2str((1:5)')); ylim(yL);
subplot(2,2,3); hold on;
	plot((1:D)-0.5,real(linreg.beta(D+(2:D+1))),'r-');
	plot(lines_x,yL'*ones(1,length(lines_x)),'-','color',0.8*[1,1,1]);
	plot((1:D)-0.5,-cx,'r.:');
	xlabel('Scale'); ylabel('Re(\beta) = cos(2\theta)'); title('Phase');
	set(gca,'xtick',[3:6:D],'xticklabel',num2str((1:5)')); ylim(yL);
subplot(2,2,4); hold on;
	plot((1:D)-0.5,imag(linreg.beta(D+(2:D+1))),'b-');
	plot(lines_x,yL'*ones(1,length(lines_x)),'-','color',0.8*[1,1,1]);
	plot((1:D)-0.5,-sx,'b.:');
	xlabel('Scale'); ylabel('Im(\beta) = sin(2\theta)'); title('Phase');
	set(gca,'xtick',[3:6:D],'xticklabel',num2str((1:5)')); ylim(yL);
graph(1); set(1,'paperposition',[0 0 20 8]); 
exportfig([figpath,'linreg_coeffs']);

return

%% histogram of true vs estimated angle
tr = linspace(-pi/2,pi/2,40);
yr = linspace(-pi/2,pi/2,40);
h = hist3([angle(y_test)/2 angle(y_linear)/2],{tr,yr});
figure(10); clf; colormap(gray(256));
	imagesc(tr,yr,h);
	
%% compute error with respect to input angle
errvec = y_linear-y_test;
cerr = real(errvec); serr = imag(errvec);
[mean(cerr) mean(serr) var(cerr) var(serr)]
[angs,inds] = sort(angle(y_test));
figure(11); clf; hold on;
	plot(angs,cerr(inds),'r-',[-pi,pi],0*[1,1],'k:'); 
	plot(angs,serr(inds)-5,'b-',[-pi,pi],-5*[1,1],'k:'); 
	plot(angs,angerr(inds)-10,'g-',[-pi,pi],-10*[1,1],'k:');
	xlim([-pi,pi]);

%% scatter plot of inputs and outputs for each method
inds = randperm(length(y_test));
N = min(2000,length(y_test));
inds = inds(1:N);
figure(12); clf; hold on; sbsz = [2,2];
	mysubplot(sbsz,1); hold on; plot(y_linear(inds),'b.'); plot(y_test(inds),'r.'); axis('equal'); title('linear');
	mysubplot(sbsz,2); hold on; plot(y_logistic(inds),'b.'); plot(y_test(inds),'r.'); axis('equal'); title('logistic');
	mysubplot(sbsz,3); hold on; plot(y_boosted(inds),'b.'); plot(y_test(inds),'r.'); axis('equal'); title('boosted');
	mysubplot(sbsz,4); hold on; plot(y_tree(inds),'b.'); plot(y_test(inds),'r.'); axis('equal'); title('forest');

return

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
