function img = gen_haar_image(inds_str)

inds = [inds_str.iimg; inds_str.iimg45_horiz; inds_str.iimg45_vert];
imsz = max(max(abs(inds(:,1:2))));

% create image buffer
img = zeros(2*imsz+1);
r = imsz+1; c = imsz+1;

inds = inds_str.iimg;
for i = 1:size(inds,1)
	img(1:r+inds(i,1),1:c+inds(i,2)) = ...
		img(1:r+inds(i,1),1:c+inds(i,2)) + inds(i,3);
end
		
inds = inds_str.iimg45_horiz;
for i = 1:size(inds,1)
	rows = r+inds(i,1); col = c+inds(i,2);
	while col>0
		img(rows,col) = img(rows,col) + inds(i,3);
		rows = max(1,rows(1)-1):min(size(img,1),rows(end)+1);
		col = col-1;
	end
end

inds = inds_str.iimg45_vert;
for i = 1:size(inds,1)
	row = r+inds(i,1); cols = c+inds(i,2);
	while row>0
		img(row,cols) = img(row,cols) + inds(i,3);
		cols = max(1,cols(1)-1):min(size(img,2),cols(end)+1);
		row = row-1;
	end
end


			