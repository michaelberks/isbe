function [data,params] = register(ref_imagefile, free_imagefile, paramfile)
%
% Pairwise registration main function
%
% Accepts filenames for reference and free images, but provides option
% to select them otherwise. Can use '' if you don't want to specify that file.
% Loads parameters from file reg_params.m if another filename is not given.

% Moves to directory of this file
myfilename = which(mfilename);
[pathname] = fileparts(myfilename);
cd(pathname);

ftype = {'*.bmp'; '*.jpg'; '*.gif'; '*.tif'; '*.png'};

% Load reference and free images and parameters
if(nargin==3)
	% Load the parameter file (an m-file)
	params = paramfile;
	% Check if image filenames are already set
	if(strcmp(free_imagefile,''))
		[filename, pathname] = uigetfile(ftype,'Free Image');
		if(isnumeric(filename))
        		break;
		end
		free_imagefile = [pathname,filename];
	end
	if(strcmp(ref_imagefile,''))
		[filename, pathname] = uigetfile(ftype,'Reference Image');
		if(isnumeric(filename))
        		break;
		end
		ref_imagefile = [pathname,filename];
	end
elseif(nargin<3)
	params = reg_params;
	if((nargin<2) | (strcmp(free_imagefile,''))) 
		[filename, pathname] = uigetfile(ftype,'Free Image');
		if(isnumeric(filename))
        		break;
		end
		free_imagefile = [pathname,filename];
	end
	if((nargin<1) | (strcmp(ref_imagefile,'')))
		[filename, pathname] = uigetfile(ftype,'Reference Image');
		if(isnumeric(filename))
        		break;
		end
		ref_imagefile = [pathname,filename];
	end
end

% If making animations, make a new directory for them
if(strcmp(params.anim,'yes'));
	% Make a new folder to put the animation figures into
	filelist = dir('./anim*');
	if(~isempty(filelist))
	    for i=1:length(filelist)
   	     f = filelist(i).name(5:7);
    	    fl(i) = str2num(f);
    	end
    	num = max(fl) + 1;
	else
	    num = 0;
	end
	name = 'anim';
	numb = num2str(num);
	if(num<10)
	    numb = strcat('00',numb);
	elseif (num<100)
	    numb = strcat('0',numb);
	end

	name = strcat(name,numb);
	params.dirname = fullfile('./',[name]);
	mkdir(params.dirname);
	params.filenum = 0;
end

% Preprocess the images to right size, intensity range, etc.
[data] = image_preprocessing(ref_imagefile,free_imagefile);
cmap = bone(256);

% Get the user to generate a region of interest on the image
[data.roi,plot_image] = region_of_interest(data.ref_image,cmap);

plot_fig = figure;, ref_axe = subplot(1,3,1);, axe = gca;, plot_surf(data.xgrid_ref,data.ygrid_ref,plot_image,cmap,axe);, set(axe,'DataAspectRatio',[1 1 1],'Xlim',[-1 1],'Ylim',[-1 1]);...
subplot(1,3,2);, axe = gca;, plot_surf(data.xgrid_free,data.ygrid_free,data.free_image,cmap,axe);, set(axe,'DataAspectRatio',[1 1 1],'Xlim',[-1 1],'Ylim',[-1 1]);

%%%%%%%% Affine registration %%%%%%%
if(params.affine==0)
	% No affine
	theta = 0;
	data.affine(1) = cos(theta);
	data.affine(2) = sin(theta);
	data.affine(3:6) = 0;
	[data.xgrid_affine,data.ygrid_affine] = affine_trans(data.xgrid_ref,data.ygrid_ref,data.affine);
else
	[data] = affine_registration(data,params);
end

% Pullback
[pullback_image] = pull_back(data.xgrid_free,data.ygrid_free,data.xgrid_affine,data.ygrid_affine,data.free_image,'median');

%%%%%%%% Non-rigid registration %%%%

%% Choose the original knot points %%
nold = 0;

if(isfield(params,'imagetype'))
	if(strcmp(params.imagetype,'skull'))
		% Pick skull points
		[xskull,yskull] = pick_skull_points(data,params.nskull,'ref');
		data.q0(1,:) = reshape(xskull,[1 params.nskull]);
		data.q0(2,:) = reshape(yskull,[1 params.nskull]);
	
		%[xskull,yskull] = pick_skull_points(data,params.nskull,'free');
		%data.q1(1,:) = reshape(xskull,[1 params.nskull]);
		%data.q1(2,:) = reshape(yskull,[1 params.nskull]);

		% Copy across and assume that affine means they are nearly in the right place!
		data.q1 = data.q0;
		if(params.perturb==1)
			data.q1 = data.q0 + 0.005*randn(size(data.q0));
		end

		nold = params.nskull;
	end
end

[disc] = discrepancy_image(data,params);

for i=1:params.nloops
	% Pick points based on discrepancy image
	[xnew,ynew,params] = pick_disc_points(data,params);
	if(~isempty(xnew))
		nnew = length(xnew);
		data.q0(1,nold+1:nold+nnew) = reshape(xnew,[1 nnew]);
		data.q0(2,nold+1:nold+nnew) = reshape(ynew,[1 nnew]);

		% Final positions of new knots using warp from old knots, otherwise MI can go up
		[xout,yout,E,data] = nrr_trans(data.q0(1,nold+1:nold+nnew),data.q0(2,nold+1:nold+nnew),data.q0(1,1:nold),data.q0(2,1:nold),data.q1(1,1:nold),data.q1(2,1:nold),params,data);
		data.q1(1,nold+1:nold+nnew) = reshape(xout,[1 nnew]);
		data.q1(2,nold+1:nold+nnew) = reshape(yout,[1 nnew]);

		% Only need to perturb if NOT warped so far!
		% Start from slight random perturbation??
		%if(params.perturb==1)
		%	data.q1(:,[nold+1:nold+nnew]) = data.q0(:,[nold+1:nold+nnew]) + 0.005*randn(size(data.q0(:,[nold+1:nold+nnew])));
		%end
	elseif(params.perturb==1)
		% Perturb anyway!		
		data.q1 = data.q0 + 0.005*randn(size(data.q0));
	end
	
	% Calculate qorig - in frame of unaffined ref image
	nknots = size(data.q0,2);
	pmap = jet(nknots);
	[qxorig,qyorig] = inverse_affine_trans(data.q0(1,:),data.q0(2,:),data.affine);
	figure(plot_fig), subplot(1,3,1);, hold on, for i=1:nknots, plot(qxorig(i),qyorig(i),'ko','MarkerFaceColor',pmap(i,:),'MarkerSize',5);, end;
	
	% Perform the non-rigid registration
	[data,params] = nrr_registration(data,params);

	% Calculate final discrepancy image
	[disc] = discrepancy_image(data,params);

	figure;, imagesc(disc);, axis image;
	nold = size(data.q0,2);
end

% Compute the pullback image
[pullback_image] = pull_back(data.xgrid_free,data.ygrid_free,data.xgrid_nrr,data.ygrid_nrr,data.free_image,'median');
figure(plot_fig);, subplot(1,3,3);, axe = gca;, plot_surf(data.xgrid_ref,data.ygrid_ref,pullback_image,cmap,axe,data.q1(1,:),data.q1(2,:),pmap);, set(axe,'DataAspectRatio',[1 1 1],'Xlim',[-1 1],'Ylim',[-1 1]);

% Output the interlaced results
interlaced1 = interlace(data.ref_image,data.free_image,params.nsquares);
interlaced2 = interlace(data.ref_image,pullback_image,params.nsquares);
figure;, subplot(1,2,1);, imagesc(interlaced1);, axis image;, subplot(1,2,2);, imagesc(interlaced2);, axis image;, colormap(bone);
