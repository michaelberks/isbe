function [haar_responses] = compute_haar_responses(img, scales)
%COMPUTE_HAAR_RESPOSES Compute responses to Haar-like filters
%   [line_strength, orientation, scale] = gaussian_2nd_derivative_line(im, scales)
%
% Inputs:
%      im - *Insert description of input variable here*
%
%      scales - *Insert description of input variable here*
%
%
% Outputs:
%      line_strength - *Insert description of input variable here*
%
%      orientation - *Insert description of input variable here*
%
%      scale - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 01-Dec-2010
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester

f_debug = false;
if (nargin==0)
	f_debug = true;
	imsz = 512; 
	img = rand(imsz,imsz);
	scales = [1 2 4 8];
% 	profile on; profile clear;
end

% pre-allocate output arguments
[r c] = size(img);
haar_responses = zeros(r, c, length(scales), 2);
haar_responses2 = haar_responses;

if f_debug
	figure(1); clf; colormap(gray(256));
end

for iscl = length(scales):-1:1
	% get parameters of first Haar-like feature
	width = 2*scales(iscl); % width is approximately 3*sigma
	n_pix00 = 4*width*width; % number of nonzero pixels
	haar_inds00 = compute_haar_inds('diag',width,width,1);

	% Parameters of second Haar-like feature (at 45deg)
	width = round(width/sqrt(2)); % reduce width for feature at 45deg
	n_pix45 = 2*4*width*width; % number of nonzero pixels
	
	% scaling to account for different numbers of pixels between filters
	norm_scl = n_pix00/n_pix45;
	haar_inds45 = compute_haar_inds('diag45',width,width,norm_scl);

	% on first pass (where filter is biggest)...
	if iscl==length(scales)
		% ...pad image with zeros...
		off = 2*width+1;
		padimg = zeros(r+2*off,c+2*off);
		padimg(off+1:end-off,off+1:end-off) = img;

		% ...and compute integral images
		intimg_str.iimg = integral_image(padimg);
		intimg_str.iimg45_horiz = integral_image_diag(padimg);
		intimg_str.iimg45_vert = integral_image_diag(padimg')';
		
		[cc,rr] = meshgrid(1+off:c+off,1+off:r+off);
		inds0	= sub2ind(size(intimg_str.iimg),rr,cc);
	end

	% get responses to Haar-like filters
	haar_responses(:,:,iscl,1) = ...
		get_haar_response(intimg_str,haar_inds00,inds0); % ~ Ixy
	haar_responses(:,:,iscl,2) = ...
		get_haar_response(intimg_str,haar_inds45,inds0); % ~ Ixx-Iyy
	
	if f_debug
		% show filters
		haar00 = gen_haar_image(haar_inds00);
		haar45 = gen_haar_image(haar_inds45);
		figure(1); colormap(gray(256));
			subplot(length(scales),2,(iscl-1)*2+1); imagesc(haar00); axis('image');
			subplot(length(scales),2,(iscl-1)*2+2); imagesc(haar45); axis('image');
			
		% check using conv2() that outputs match
		output00 = conv2(img,haar00,'same');
		err00 = norm(output00-haar_responses(:,:,iscl,1),'fro');
		output45 = conv2(img,haar45,'same');
		err45 = norm(output45-haar_responses(:,:,iscl,2),'fro');
		if (err00>1e-6 | err45>1e-6)
			error('Outputs do not match');
		end
		
		% a few more checks
		if ( abs(sum(abs(haar00(:)))-sum(abs(haar45(:))))>1e-6 )
			error('Filters do not have correct weights');
		end
	end
end

if f_debug
% 	profile off; profile report;
	clear;
end


function resp = get_haar_response(intimg_str,haar_inds,inds0)
% get response to Haar-like filter

% indices and their offsets for row and column
roff	= 1; 
coff	= size(intimg_str.iimg,1);

resp	= zeros(size(inds0));
if ~isempty(haar_inds.iimg)
	indoffs	=	roff*haar_inds.iimg(:,1) + ...
				coff*haar_inds.iimg(:,2);
	for i = 1:size(haar_inds.iimg,1)
		inds = inds0+indoffs(i);
		resp = resp + intimg_str.iimg(inds)*haar_inds.iimg(i,3);
	end
end

if ~isempty(haar_inds.iimg45_horiz)
	indoffs	=	roff*haar_inds.iimg45_horiz(:,1) + ...
				coff*haar_inds.iimg45_horiz(:,2);
	for i = 1:size(haar_inds.iimg45_horiz,1)
		inds = inds0+indoffs(i);
		resp = resp + intimg_str.iimg45_horiz(inds)*haar_inds.iimg45_horiz(i,3);
	end
end

if ~isempty(haar_inds.iimg45_vert)
	indoffs	=	roff*haar_inds.iimg45_vert(:,1) + ...
				coff*haar_inds.iimg45_vert(:,2);
	for i = 1:size(haar_inds.iimg45_vert,1)
		inds = inds0+indoffs(i);
		resp = resp + intimg_str.iimg45_vert(inds)*haar_inds.iimg45_vert(i,3);
	end
end
