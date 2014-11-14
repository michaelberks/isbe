% clear;

% load distances
% outpath = 'U:\projects\mammography\data\histograms\';
load(fullfile(outpath,'vq_dists_normals.mat'));
load(fullfile(outpath,'vq_dists_abnormals.mat'));


% get rid of those examples used for training
vq_dists_normals(isnan(vq_dists_normals)) = [];

n_normal = length(vq_dists_normals);
n_abnormal = length(vq_dists_abnormals);

dh = 0.1; d_range = 0:dh:1;
h_normals = hist(vq_dists_normals,d_range);
h_abnormals = hist(vq_dists_abnormals,d_range);

display([mean(vq_dists_normals) mean(vq_dists_abnormals)]);
legstrs = { sprintf('Normals (n=%i, \\mu=%0.2f)',n_normal,mean(vq_dists_normals)),...
			sprintf('Abnormals (n=%i, \\mu=%0.2f)',n_abnormal,mean(vq_dists_abnormals)) };

figure(1); clf; hold on;
% normal samples
subplot(2,2,1);
	plot(vq_dists_normals([1,1],:),[zeros(1,n_normal); ones(1,n_normal)],'b-');
	axis([0,1,0,1.2]);
	xlabel('L-R distance'); 
	title('Normal samples');
% abnormal samples
subplot(2,2,3);
	plot(vq_dists_abnormals([1,1],:),[zeros(1,n_abnormal); ones(1,n_abnormal)],'r-');
	axis([0,1,0,1.2])
	xlabel('L-R distance'); 
	title('Abnormal samples');
% histograms
subplot(2,2,2); box on;
	plot(d_range,h_normals/(n_normal*dh),'b-',...
		 d_range,h_abnormals/(n_abnormal*dh),'r-');
	xlabel('L-R distance'); 
	title('Histogram');
	legend(legstrs);
% cumulative histograms
subplot(2,2,4); box on;
	plot([0 sort(vq_dists_normals)],(0:n_normal)/n_normal,'b-',...
		 [0 sort(vq_dists_abnormals)],(0:n_abnormal)/n_abnormal,'r-');
	axis([0,1,0,1])
	xlabel('L-R distance');
	title('Cumulative Histogram');
	legend(legstrs,'location','southeast');

% outpath = 'U:\projects\mammography\figs\vq_comparison\';
exportfig(fullfile(outpath,'graphs'));
