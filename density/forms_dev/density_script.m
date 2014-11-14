% Load in the components of the form

marker1 = imread('K:\isbe\density\form\marker1.bmp');
marker2 = imread('K:\isbe\density\form\marker2.bmp');
marker3 = imread('K:\isbe\density\form\marker3.bmp');

%line_scale_im = imread('K:\isbe\density\form\line_scale.bmp');
%line_scale = line_scale_im > 0;

line_scale = true(125, 600);
line_scale(62:64, 101:500) = 0;

patient_id = double(rgb2gray(imread('k:\isbe\density\form\patient_id.bmp')));
%
save K:\isbe\density\form\components marker* line_scale 
%

% Put the componenets together to build the complete form
form = [marker1 ones(50, 550);...
        ones(100, 600);
        line_scale; line_scale; line_scale; line_scale;
        ones(50, 600);
        marker2 ones(50, 500) marker3];
imwrite(form, 'K:\isbe\density\form\form_image.bmp');
%%

%Try and find marker and marker90

C1 = normxcorr2(marker1, form);
C2 = normxcorr2(marker2, form);
C3 = normxcorr2(marker3, form);

[dummy, ind1] = max(C1(:));
[pos1(1) pos1(2)] = ind2sub(size(C1), ind1);
pos1 = pos1 - [25 25];

[dummy, ind2] = max(C2(:));
[pos2(1) pos2(2)] = ind2sub(size(C2), ind2);
pos2 = pos2 - [25 25];

[dummy, ind3] = max(C3(:));
[pos3(1) pos3(2)] = ind2sub(size(C3), ind3);
pos3 = pos3 - [25 25];
%
figure; imagesc(form); axis image; colormap(gray(256)); hold on;

plot(pos1(2), pos1(1), 'rx', 'MarkerSize', 10);
plot(pos2(2), pos2(1), 'bx', 'MarkerSize', 10);
plot(pos3(2), pos3(1), 'gx', 'MarkerSize', 10);

%%
% Marvellous, now find the start and end points of the line scales

line_start1 = pos0 + [110 35];
line_end1 = line_start1 + [0 400];
plot([line_start1(2) line_end1(2)], [line_start1(1) line_end1(1)], 'g:');

line_start2 = pos0 + [210 35];
line_end2 = line_start2 + [0 400];
plot([line_start2(2) line_end2(2)], [line_start2(1) line_end2(1)], 'g:');

line_start3 = pos0 + [310 35];
line_end3 = line_start3 + [0 400];
plot([line_start3(2) line_end3(2)], [line_start3(1) line_end3(1)], 'g:');

line_start4 = pos0 + [410 35];
line_end4 = line_start4 + [0 400];
plot([line_start4(2) line_end4(2)], [line_start4(1) line_end4(1)], 'g:');

%%
patient_id = double(rgb2gray(imread('k:\isbe\density\form\patient_id.bmp')));
%
patient_id_orig = cell(8, 10);

[xx yy] = meshgrid(5:20, 3:18);
%
figure; image(patient_id); axis image; colormap(gray(256)); hold on;
for row = 1:8
    for col = 1:10
        patient_id_orig{row,col} = [xx(:)+(col-1)*23 yy(:)+(row-1)*19];
        plot(patient_id_orig{row,col}(:,1),  patient_id_orig{row,col}(:,2), 'r.');
        
        digit_pix = interp2(patient_id, patient_id_orig{row,col}(:,1), patient_id_orig{row,col}(:,2));
        display(sum(digit_pix));
    end
end
%%
for ii = 1:8
    figure; 
    image(reshape(interp2(patient_id, patient_id_orig{ii,ii}(:,1), patient_id_orig{ii,ii}(:,2)), 16, 16)); 
    axis image; colormap(gray(256));
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute binary number from base 10 into a 32 digit binary number and then
% generate barcode
pma_orig = round((2^31-1)*rand);
pma = pma_orig;
pma_bin = zeros(1,32);
for ii = 31:-1:0
    pma_bin(32 - ii) = 2^ii <= pma;
    pma = rem(pma, 2^ii);
end
% pma_barcode = kron(pma_bin, ones(16, 4));
% figure; imagesc(pma_barcode); axis image; colormap(gray(256));
% imwrite(pma_barcode, 'K:\isbe\density\form\barcode.bmp');

%reconstruct pma
pma_recon = 0;
for ii = 31:-1:0
    pma_recon = pma_recon + (pma_bin(32 - ii)* 2^ii);
end

pma_orig - pma_recon
pma_orig