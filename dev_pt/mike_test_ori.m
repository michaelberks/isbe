clc;

imtype = 'real';

nty = '';
% nty = '_90';

contrast = 4; contstr = '';

%% Make orientation images:
outpath = [asymmetryroot,'data\synthetic_lines\ori_',imtype,'512',nty,'\labels'];
if ~exist(outpath,'dir'), mkdir(outpath); end

disp('Generating images...');
for ii = [1 31 61]
	inpath	= sprintf('%sdata/synthetic_backgrounds/%s512/test/bg%05i.mat',...
										asymmetryroot,imtype,ii);
	bg = u_load(inpath);

	for jj = 1:6:180
		outpath	= sprintf('%sdata/synthetic_lines/ori_%s512%s/image_%03i_%03i%s.mat',...
											asymmetryroot,imtype,nty,jj,ii,contstr);
		lblpath	= sprintf('%sdata/synthetic_lines/ori_%s512%s/labels/label_%03i_%03i%s.mat',...
											asymmetryroot,imtype,nty,jj,ii,contstr);

		% if images exist then skip
		if (exist(outpath,'file') && exist(lblpath,'file')), continue; end

		[bar_image, label, label_centre, label_orientation] = ...
			create_sin_bar(3, contrast, jj, 512, 512, 0.5, 256, 256);

		if strcmp(imtype,'real')
			if strcmp(nty,'_90')
				bar_real_90 = bar_image + rot90(bg);
				save(outpath,'bar_real_90');
			else
        bar_real = bar_image + bg;
				save(outpath,'bar_real');
			end
		else
			bar_syn = bar_image + bg;
			save(outpath,'bar_syn');
		end
		
		save(lblpath,'label','label_centre','label_orientation');
	end
end

%% Compute orientation maps for g2 method (RF computed on hydra)
outpath = [asymmetryroot,'data\synthetic_lines\ori_',imtype,'512',nty,'\results\g2_ori'];
if ~exist(outpath,'dir'), mkdir(outpath); end

disp('Computing G2D orientations...');
for jj = 1:6:180
	for ii = [1 31 61]
		outpath	= [	asymmetryroot,'data\synthetic_lines\ori_' imtype '512',nty,'\results\g2_ori\image'...
								zerostr(jj,3) '_' zerostr(ii,3),contstr,'_ori.mat'];

		% if orientation already computed then skip
		if exist(outpath,'file'), continue; end
		
		inpath	= [	asymmetryroot,'data\synthetic_lines\ori_',imtype,'512',nty,'\image_'...
								zerostr(jj,3) '_' zerostr(ii,3),contstr,'.mat'];
		test_image = load_uint8(inpath);

		[orientation_map] = karssemeijer_line_detection(...
				test_image,...
				'line_scales', [1 2 4 8],...
				'grad_scale', 10,...
				'grad_ori_thresh', pi/6,...
				'grad_strength_thresh', 25,...
				'line_strength_thresh', 0,...
				'binary_map', 1);

		save_uint8(outpath, orientation_map);
	end
end

%% Copy orientation maps for RF method from hydra
rfid = '238470'; % dtcwt RF
% rfid = '233142'; % g2d RF
disp('Copying RF orientations...');
inpath	= ['A:\data\synthetic_lines\ori_',imtype,'512',nty,'\results\',rfid]; 
outpath = [asymmetryroot,'data\synthetic_lines\ori_',imtype,'512',nty,'\results\',rfid]; 
if ~exist(outpath,'dir'), copyfile(inpath,outpath); end
	
%% Look at properties of the RF used to predict orientations
load([asymmetryroot,'data/line_orientation_rfs/',rfid,'/random_forest.mat']);
rf = random_forest;
rf.tree_root = strrep(rf.tree_root,...
											'/san/images/asym/',asymmetryroot);
[rf_hist_in,rf_hist_out] = rf_training_hist(rf);
load([asymmetryroot,'data/line_orientation_rfs/',rfid,'/sampling_args.mat']);

%% Produce error plots of the various image for the g2/RF methods
% Orientations
xy = repmat(-256:255, 512, 1);
circle_mask = xy.^2 + xy'.^2 < 128^2;
mean_errors = [];

confmat_rf = zeros(180,180); % confusion matrix
confmat_g2 = zeros(180,180); % confusion matrix

disp('Computing errors...');
trng = 1:6:180;
for ori_gt = trng
    ori_errors = [];

		for ii = [1 31 61]
        ori_map_rf = load_uint8([asymmetryroot,'data\synthetic_lines\ori_' imtype '512',nty,'\results\',rfid,'\image'...
            zerostr(ori_gt,3) '_' zerostr(ii,3) '_class.mat']);
        ori_map_g2 = load_uint8([asymmetryroot,'data\synthetic_lines\ori_' imtype '512',nty,'\results\g2_ori\image'...
            zerostr(ori_gt,3) '_' zerostr(ii,3) '_ori.mat']);
        label = load([asymmetryroot,'data\synthetic_lines\ori_' imtype '512',nty,'\labels\label_'...
            zerostr(ori_gt,3) '_' zerostr(ii,3) '.mat']);

        label_centre = label.label_centre & circle_mask;

				% saved outputs do not have double angles
				ori_map_rf2 = abs(ori_map_rf).*exp(complex(0,2*angle(ori_map_rf)));
				
				% compute orientation errors (ori_gt is in degrees)
				ori_errors = [ori_errors; ...
					ori_error(ori_map_rf2(label_centre),ori_gt,'correct','deg') ...
					ori_error(ori_map_g2(label_centre),ori_gt,'correct','deg')]; %#ok

				% compute histogram over 
				hst = histc(angle(ori_map_rf(label_centre))*(180/pi),-90.5:90.5)';
					hst(1) = hst(1)+hst(end-1); hst = hst(1:end-2);
					confmat_rf(ori_gt,[90:180,1:89]) = confmat_rf(ori_gt,[90:180,1:89]) + hst;

				hst = histc(ori_map_g2(label_centre),-45.5:135.5)';
					hst(1) = hst(1)+hst(end-1); hst = hst(1:end-2);
					confmat_g2(ori_gt,[135:180,1:134]) = confmat_g2(ori_gt,[135:180,1:134]) + hst;
		end
% 		confmat_rf(ori_gt,:) = confmat_rf(ori_gt,:)/sum(confmat_rf(ori_gt,:));
% 		confmat_g2(ori_gt,:) = confmat_g2(ori_gt,:)/sum(confmat_g2(ori_gt,:));
		
    mean_errors(end+1,:) = mean(ori_errors*180/pi); %#ok
end

figure(1); clf; hold on;
	plot(trng, mean_errors(:,1), 'g', 'linewidth', 2);
	plot(trng, mean_errors(:,2), 'c', 'linewidth', 2);
	axis([-4 180 0 60]);
	title('Orientation errors for lines of varying orientation');
	ylabel('Mean absolute orientation error (degrees)');
	xlabel('Line orientation (degrees)');
	set(gca, 'xtick', trng, 'xticklabel', trng);

confmat_rf = confmat_rf(1:6:180,:);
confmat_g2 = confmat_g2(1:6:180,:);

figure(2); clf;
	sbsz = [2,2];
	mysubplot(sbsz);
		imagesc(1:180,1:6:180,confmat_g2); colormap(gray(256)); axis('square','xy');
% 		hold on; plot(confmat_g2*(1:180)',1:6:180,'r-');
		xlabel('Estimated \theta'); ylabel('True \theta');
	mysubplot(sbsz);
		imagesc(1:180,1:6:180,confmat_rf); colormap(gray(256)); axis('square','xy');
% 		hold on; plot(confmat_rf*(1:180)',1:6:180,'r-');
		xlabel('Estimated \theta'); ylabel('True \theta');

% figure(3); clf;
% 	sbsz = [1,2];
	mysubplot(sbsz); 
		plot(sum(confmat_g2)/sum(confmat_g2(:))); 
		axis([0,180,0,0.015]);
		xlabel('Estimated \theta'); ylabel('Prob. dens.');
	mysubplot(sbsz); hold on;
		plot(rf_hist_out,'r-');
		plot(sum(confmat_rf)/sum(confmat_rf(:)),'b-');
		axis([0,180,0,0.015]);
		xlabel('Estimated \theta'); ylabel('Prob. dens.');
 		legend({'Potential output','Actual output'},'location','northwest');
		set(gca,'box','on');
	
outroot = 'U:\matlab\figs\mammography\';
exportfig([outroot,'confmats_',rfid,nty,contstr]);
