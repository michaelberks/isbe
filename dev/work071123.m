list = dir('C:\isbe\dev\background\images\mass\', 'mass*');
%%
new_patches = [];
idx = 1:53;
for jj = 1:20
    mass = uint8(new_masses(jj).mass_ROI);
    
    for ii = 1:length(idx);
        patch = imread(['C:\isbe\dev\background\images\mass\', list(idx(ii)).name]);
        if all(size(mass) <= size(patch))
            [rows_m cols_m] = size(mass);
            [rows_p cols_p] = size(patch);
            
            r_start = 1+ round((rows_p - rows_m)/2);
            c_start = 1+ round((cols_p - cols_m)/2);
            r_end = r_start + rows_m - 1;
            c_end = c_start + cols_m - 1;
            
            patch(r_start:r_end, c_start:c_end) = ...
                patch(r_start:r_end, c_start:c_end) + mass;
            figure; imagesc(patch); colormap(gray(256)); axis image;
            new_patches(end+1).roi = patch;
            new_patches(end).idx = idx(ii);
            clear mass; clear patch;
            idx(ii) = [];
            break;
        end
        clear patch;
    end
    clear mass;
end
%%
load C:\isbe\dev\files\bg_files.mat
for ii = 1:10
   figure; 
   subplot(1,2,1); image(new_patches(ii).roi); colormap(gray(256)); axis image;
   load(['C:\isbe\dev\masses\', bg_files(new_patches(ii).idx).name]);
   subplot(1,2,2); image(mass.mass_ROI); colormap(gray(256)); axis image;
end
%%
load C:\isbe\dev\files\u_files.mat
load C:\isbe\dev\models\model_w500_50K.mat
%%
model_errors_loo2(mass_model, u_files1, 'save_recon', 'C:\isbe\dev\recon\');
%%
for ii = 1:20
    load(['C:\isbe\dev\masses\', u_files1(ii).name]);
    r1 = load(['C:\isbe\dev\recon\recon_', zerostr(ii,3)]);
    figure; 
    subplot(1,3,1); imagesc(mass.subtract_ROI); colormap(gray(256)); axis image;
	subplot(1,3,2); imagesc(r2.old_shape_ROI); colormap(gray(256)); axis image;
	subplot(1,3,3); imagesc(r2.new_shape_ROI); colormap(gray(256)); axis image;
end