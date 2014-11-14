function test_nrr_warps

% Grid of points
inputqx = [];
inputqy = [];
load 'E:\Diffeoflow\big_grid_points.mat'

xin = inputqx;
yin = inputqy;
clear inputqx inputqy;

[qx0,qy0] = meshgrid([-0.5 0 0.5],[-0.5 0 0.5]);

qx0 = [qx0(:); -0.9; 0.9; 0; 0];
qy0 = [qy0(:);  0; 0; 0.9; -0.9];

nknots = length(qx0(:))

nset = 5;

amp = 0.3;

types = green;
ntypes = length(types);

for i=1:nset
	theta = 2*pi*rand(size(qx0));
	r = amp*rand(size(qx0));

	qx1 = qx0 + r.*cos(theta);
	qy1 = qy0 + r.*sin(theta);

	% Add fixed points at corners
	qx1(nknots-4:nknots) = qx0(nknots-4:nknots);
	qy1(nknots-4:nknots) = qy0(nknots-4:nknots);



	for j=1:ntypes
		params.green = char(types(j));

		[xout,yout,E,data] = nrr_trans(xin,yin,qx0,qy0,qx1,qy1,params,data);

		figure;, subplot(1,2,1), plot(xin,yin,'k.','MarkerSize',1);, hold on, plot(qx0,qy0,'ro','MarkerFaceColor','r','MarkerSize',5);,...
		plot(qx1,qy1,'go','MarkerFaceColor','g','MarkerSize',3);, set(gca,'DataAspectRatio',[1 1 1]);
		subplot(1,2,2), plot(xout,yout,'k.','MarkerSize',1);, hold on, plot(qx0,qy0,'ro','MarkerFaceColor','r','MarkerSize',5);,...
		plot(qx1,qy1,'go','MarkerFaceColor','g','MarkerSize',3);, xlabel(char(params.green));,, set(gca,'DataAspectRatio',[1 1 1]);
	end
end	
