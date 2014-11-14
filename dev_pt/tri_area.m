imgpath = [asymmetryroot,'data\mammograms\2004_screening_processed\mass_roi\png\'];
imgdir	= dir([imgpath,'*.png']);

oripath	= 'A:\data\orientation_maps\rf\2004_screening_processed\mass_roi\';
oridir	= dir([oripath,'*_roi.mat']);

for i = 1:1
	if ~exist('img','var'),
		img = imread([imgpath,imgdir(i).name]);
		ori = load_uint8([oripath,oridir(i).name],true);
	end

	% show the mammogram ROI
	figure(1); clf; colormap(gray(256));
	subplot(2,2,1);
		image(uint8(img));
		axis('off','image');
		
	% compute a colour image from the orientation map and show it
	ori_hsv(:,:,1) = (angle(ori)/pi) + 0.5;
	ori_hsv(:,:,2) = 1;
	ori_hsv(:,:,3) = abs(ori);
	ori_hsv = hsv2rgb(ori_hsv);
	subplot(2,2,2);
		image(ori_hsv);
		axis('off','image');
	
	% compute the distance of all pixels from the mass centre
	mass_c2 = [320,378];
	[xx,yy] = meshgrid( (1:size(ori,2))-mass_c2(1),...
						(1:size(ori,1))-mass_c2(2) );
	rad = sqrt(xx.*xx+yy.*yy);
	
	% define thresholds for dispersion/strength
	thr_list = 0:0.1:0.9;
% 	thr_list = 0;
	
	% define sigmas for synthetic experiments
% 	sig_list = -1; % use estimated directions
	sig_list = [-1 logspace(0,-2,9)]; % add noise to ideal directions

	ntrials = 50;

	% preallocate space
	stats = zeros(length(thr_list),length(sig_list),ntrials);
	
	for ithr = 1:length(thr_list)
		% find all pixels within the donut that are above the threshold
		thr = thr_list(ithr);
		mask = rad>100 & rad<200 & abs(ori)>thr;
		mask_inds = find(mask)';
		
		% convert to row/column indices
		[xx,yy] = meshgrid( 1:size(ori,2),1:size(ori,1) );
		cc = xx(mask_inds); rr = yy(mask_inds);

		for isig = 1:length(sig_list)
			sig = sig_list(isig);
			
			if (ithr>1 && isig>1), break; end
			
			% compute the histogram ntrials times
			for trial = 1:ntrials
				log_areas = nan(1,1000);

				for p = 1:length(log_areas)
					% randomly sample 3 pixels from the valid set

					% randperm is ok for small pixel sets...
			% 		inds = randperm(numel(ori_tmp)); inds = inds(1:3);
			% 		inds = randperm(length(mask_inds)); inds = inds(1:3);

					% ...but this is faster for bigger sets
					inds = ceil(rand(1,3)*length(mask_inds));
					while length(unique(inds))<3
						inds = ceil(rand(1,3)*length(mask_inds));
					end

					% if orientation is not defined for any of the three points
					% then ignore this sample (should really resample)
					if any(ori(mask_inds(inds))==0), continue; end

					% locations of sampled points
					v = [cc(inds); rr(inds)];
					
					% direction vectors
					if sig==-1
						% real
						d = ori(mask_inds(inds))./abs(ori(mask_inds(inds))); 
						d = [real(d); -imag(d)]; % y-flipped due to coord frame
					else
						% synthetic (ideal direction + noise)
						t = atan2(v(2,:)-mass_c2(2),v(1,:)-mass_c2(1)) + sig*randn(1,3); 
						d = [cos(t); sin(t)];
					end
					
					% compute corners of triangle
					c = zeros(2,3);	i2 = [1,2,3,1];
					for i = 1:3
						ns = null([v(:,i2(i))-v(:,i2(i+1)) d(:,i2(i)) d(:,i2(i+1))]);
						c(:,i) = v(:,i)+ns(2)/ns(1)*d(:,i);
					end

					% compute area of resulting triangle
					c2 = c(:,[2,3]) - c(:,[1,1]);
					log_areas(p) = log(0.5*sqrt(det(c2'*c2)));
					if isinf(log_areas(p)) | ~isreal(log_areas(p))
						log_areas(p) = NaN; 
					end
				end % for p
				
				% lose NaNs and compute log-area
				log_areas(isnan(log_areas)) = [];
				stats(ithr,isig,trial) = mean(log_areas);
			end % for trial
		end % for isig
	end % for ithr
	
	subplot(2,2,3); cla; hold on;
	set(gca,'color',[0,0,0]);
	for ii = 1:length(inds)
		i = inds(ii);
% 		p1 = complex(cc(i),rr(i)) + 1000*ori(mask_inds(i))/abs(ori(mask_inds(i)));
% 		p2 = complex(cc(i),rr(i)) - 1000*ori(mask_inds(i))/abs(ori(mask_inds(i)));
		p1 = complex(cc(i),rr(i)) + 1000*complex(d(1,ii),d(2,ii));
		p2 = complex(cc(i),rr(i)) - 1000*complex(d(1,ii),d(2,ii));
		plot(real([p1; p2]),imag([p1; p2]),'-',...
				'color','g');%abs(ori(mask_inds(i)))*[0 1 0]);
	end
	rectangle('position',[mass_c2-200,400*[1,1]],'curvature',[1,1],'edgecolor','y');
	rectangle('position',[mass_c2-100,200*[1,1]],'curvature',[1,1],'edgecolor','y');
	plot(c(1,:),c(2,:),'r.');
	plot(cc(inds),rr(inds),'c.');
	axis('equal','ij',[0,size(ori,2),0,size(ori,1)]);
	
	subplot(2,2,4); 
		hist(log_areas,0:20);
		axis([-1,21,0,length(log_areas)/4]);
	display([mean(log_areas)])
		
	subplot(2,2,1);
		rectangle('position',[mass_c2-200,400*[1,1]],'curvature',[1,1],'edgecolor','y');
		rectangle('position',[mass_c2-100,200*[1,1]],'curvature',[1,1],'edgecolor','y');
		
	figure(2);
	subplot(2,1,1);
		stats1 = reshape(stats(1,:,:),[length(sig_list),ntrials])';
		boxplot(stats1);
		set(gca,'xticklabel',num2str(sig_list','%5.2f'));
		xlabel('Noise Sigma (rad)');
	subplot(2,1,2);
		stats1 = reshape(stats(:,1,:),[length(thr_list),ntrials])';
		boxplot(stats1);
		set(gca,'xticklabel',num2str(thr_list'));
		xlabel('Threshold');
		
	figpath = 'U:\projects\mammography\figs\kmaps\';
	set(gcf,'paperposition',[0,0,16,20]);
	exportfig([figpath,'tri_area']);
end

