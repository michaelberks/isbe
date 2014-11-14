% test to investigate the sensitivity of DTCWT coefficients to noise
% if we can estimate the uncertainty in DTCWT coefficients then we can
% sample probabilistically from a tree

fhw = 4; % filter half-width

inpath = 'U:\projects\mammography\data\synthetic_lines2\';

%% choose points to sample
load([inpath,'labels/label001.mat']);

% mask out pixels near boundary
label([1:fhw,end-fhw:end],:) = false;
label(:,[1:fhw,end-fhw:end]) = false;

%% generate noisy versions of input image
load([inpath,'image001.mat']);

% smooth the image first
x = linspace(-3,3,2*fhw+1);
g = exp(-0.5*x.*x); g = g/sum(g);
g2 = g'*g;
smooth_image = conv2(g,g,test_image,'same');

% get difference image
diff_image = test_image-smooth_image;
diff_image = diff_image(fhw+1:end-fhw,fhw+1:end-fhw);

% compute noise statistics
if 0
	base_image = smooth_image;
	mn = mean(diff_image(:));
	sd = std(diff_image(:));
else
	base_image = test_image;
	mn = 0;
	sd = sqrt(base_image)/10;
end

% set sampling arguments
sampling_args = struct(	'win_size',1,...
						'feature_type','conj');
					
% clear coeffs;
figure(1); clf; colormap(gray(256));
if ~exist('coeffs','var')
	% choose N pixels
	N = 100;
	inds = find(label);
	inds = inds(randperm(length(inds)));
	[rows cols] = ind2sub(size(label),inds(1:N));

	% get initial dtcwt coefficients
	dt = compute_dual_tree(test_image,5,0);
	coeffs0 = sample_dt_data(dt,rows,cols,sampling_args);
	
	coeffs = [];
	for i = 1:500
		% generate a noise image
		noise_image = mn + sd.*randn(size(base_image));
		new_image = base_image + noise_image;
% 		imagesc(new_image(fhw+1:end-fhw,fhw+1:end-fhw)); axis('image','off');

		% get dtcwt coefficients
		dt = compute_dual_tree(new_image,5,0);
		coeffs(:,:,i) = sample_dt_data(dt,rows,cols,sampling_args);
	end
end

ori_inds = [12,8,4,1,5,9, 11,7,3,2,6,10];
th = [15,45,75,105,135,165]*pi/180;
star_x = 32*[cos(th); -cos(th)];
star_y = 32*[sin(th); -sin(th)];
for level = 2
	for point = 1:size(coeffs,1)
		figure(1); clf; hold on;
			imagesc(new_image);
			plot(cols(point),rows(point),'r.');
			plot(cols(point)+star_x,rows(point)+star_y,'r-');
			axis('ij','image','off');

		figure(2); clf;
		for half = 0:1
			params = (30*half)+(6*level)+(1:6);
			switch half,
				case 0,	
					mx = max(max(coeffs(point,params,:)));
					orivar = var(coeffs0(point,params))
				case 1,	
					mx = pi;
			end
			
			for ori = 1:6
				param = params(ori);
				c_off = coeffs0(point,param);
				
				switch half,
					case 0,	rng = c_off+[-2,2];
					case 1,	rng = [0,mx];
				end
			
				bc = linspace(rng(1),rng(2),51);
				h = hist(coeffs(point,param,:),bc);
				
				mysubplot([3,4],ori_inds((6*half)+ori)); hold on;
					plot(bc,h,'b-');
					plot([c_off,c_off],[0,100],'r-');
				axis([0,mx,0,100]);
				title(sprintf('Level %i, Ori %i, Point %i',level+1,ori,point));
			end
		end
		
		pause;
	end
end
return

