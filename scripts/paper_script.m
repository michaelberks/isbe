%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Paper script this script uitilises parts of other module scripts to
% generate figure and results that are useful for paper write-ups, formal
% reports etc.

%%
% Analysing regenerated LOO regions

load C:\isbe\dev\files\u_files.mat
%%
for ii = 81:101
    load(['C:\isbe\dev\masses\', u_files1(ii).name]);
    load(['C:\isbe\dev\recon\recon_', zerostr(ii, 3)]);
    figure;
    subplot(1,3,1);
    imagesc(mass.subtract_ROI); colormap(gray(256)); axis image;
    xlabel('Original');
    subplot(1,3,2);
    imagesc(old_shape_ROI); colormap(gray(256)); axis image;
    xlabel('Individual models');
    subplot(1,3,3);
    imagesc(new_shape_ROI); colormap(gray(256)); axis image;
    xlabel('Combined models');
    clear mass *_shape_ROI
end
%%
% the nice ones are
nice_ones = [3 21 25 45 50 63 72 80 89 92 99];
for ii = 1:length(nice_ones)
    kk = nice_ones(ii);
    load(['C:\isbe\dev\masses\', u_files1(kk).name]);
    load(['C:\isbe\dev\recon\recon_', zerostr(kk, 3)]);
    figure;
    subplot(1,3,1);
    imagesc(mass.subtract_ROI); colormap(gray(256)); axis image;
    xlabel('Original');
    subplot(1,3,2);
    imagesc(old_shape_ROI); colormap(gray(256)); axis image;
    xlabel('Individual models');
    subplot(1,3,3);
    imagesc(new_shape_ROI); colormap(gray(256)); axis image;
    xlabel('Combined models');
    clear mass *_shape_ROI
end
%%
% my favourite = 45;
load(['C:\isbe\dev\masses\', u_files1(45).name]);
load('C:\isbe\dev\recon\recon_045');

figure;
subplot(2,2,1);
imagesc(mass.subtract_ROI); colormap(gray(256)); axis image;
xlabel('Original');
subplot(2,2,2);
imagesc(new_shape_ROI); colormap(gray(256)); axis image;
xlabel('Individual models');
subplot(2,2,3);
imagesc(new_shape_ROI+double(mass.mass_ROI) - mass.subtract_ROI); colormap(gray(256)); axis image;
xlabel('Combined models');
subplot(2,2,4);
imagesc(double(mass.mass_ROI)); colormap(gray(256)); axis image;
xlabel('Combined models');





    