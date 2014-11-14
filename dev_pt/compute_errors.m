function compute_errors

clc; clear all;

estroot	= ['A:\data\synthetic_lines\real512\results\'];
rootdir = dir(estroot);

fid = fopen([asymmetryroot,'data\synthetic_lines\real512\results\summary.xls'],'w');
for i_rf = 1:length(rootdir)
	rf = str2num(rootdir(i_rf).name);
	
	% ignore anything that's not a valid number
	if isempty(rf), continue; end
	rf = num2str(rf);
	
	% process all subdirs (including '.')
	subdirs = dir([estroot,rf]);
	for i_subdir = 1:length(subdirs)
		es = process_dir([estroot,rf,'\',subdirs(i_subdir).name]);
		if ~isempty(es)
			summvec = [	es.mean,es.sd,...
						es.abs_mean,es.abs_median,es.abs_range,...
						es.abs_percentiles([10,25,75,90]) ];
			fprintf(fid,'%s/%s\t',rf,subdirs(i_subdir).name);
			fprintf(fid,'%f\t',summvec);
			fprintf(fid,'\n');
		end
	end
end
fclose(fid);
fclose all;


function err_stats = process_dir(dirname)

err_stats = [];
lblpath	= 'A:\data\synthetic_lines\real512\labels';
files	= dir([dirname,'\*_class.mat']);
n_img	= length(files); 
if (n_img==0), return; end

% n_img = 1;
errors	= cell(1,n_img);
for img = 1:n_img
	lbl = load(sprintf('%s\\label%03d.mat',lblpath,img));
	est = load_uint8(sprintf('%s\\image%03d_class.mat',dirname,img));

	inds = (lbl.label==1) & (lbl.label_centre);
	lbl = lbl.label_orientation(inds);
	est = est(inds);

	% RF regression function gives true angle so we need to double it
	est = abs(est).*exp(sqrt(-1)*2*angle(est));

	errors{img} = ori_error(lbl,est,'correct','deg') * 180/pi;
end

% stats of angular error
errvec = cat(1,errors{:});
err_stats.mean = mean(errvec);
err_stats.sd = std(errvec);

% stats of absolute angular error
errvec = sort(abs(errvec));
err_stats.abs_mean = mean(errvec);
err_stats.abs_median = median(errvec);
err_stats.abs_range = errvec([1,end])';
pcs = round([0.01:0.01:1]*length(errvec));
err_stats.abs_percentiles = errvec(pcs)';

% For the first task of job, dump the username, forest arguments and
% sampling arguments to a text file.
dirname = strrep(dirname,'A:\',asymmetryroot);
if ~exist(dirname,'dir'), mkdir(dirname); end

% copyfile([oripath,'\args.txt'],estpath);
filename = [dirname,'/err_stats.txt'];
fid = fopen(filename,'w');
	fprintf(fid,'User = %s\n',get_username());
	fprintf(fid,'%s\n',evalc('err_stats'));
	fprintf(fid,'%s\n',evalc('percentiles = [10,25,50,75,90; err_stats.abs_percentiles([10,25,50,75,90])]'''));
fclose(fid);
