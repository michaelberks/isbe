%load in markers
marker1 = imread('C:\isbe\density\form\components\marker1.bmp');
marker2 = imread('C:\isbe\density\form\components\marker2.bmp');
marker3 = imread('C:\isbe\density\form\components\marker3.bmp');

% Build lines
line_scale = true(125, 600);
line_scale(62:64, 101:500) = 0;

%Load in patient ID box
patient_id_box = double(rgb2gray(imread('C:\isbe\density\form\components\patient_id_box.bmp')));
patient_id_box = patient_id_box / 255;

%Make a random barcode (perhaps in future we'll make this script a method
%with the pma as an argument). For now lets just randomly generate one
pma_orig = round((2^31-1)*rand);
pma = pma_orig;
pma_bin = zeros(1,32);
for ii = 31:-1:0
    pma_bin(32 - ii) = 2^ii <= pma;
    pma = rem(pma, 2^ii);
end
pma_barcode = kron(pma_bin, ones(20, 4));

% Save the components
save C:\isbe\density\form\components\components marker* line_scale patient_id_box 

% Put the componenets together to build the complete form
form = [marker1 ones(50, 550);...
        ones(200, 600);
        line_scale; line_scale; line_scale; line_scale;
        ones(50, 600);
        marker2 ones(50, 500) marker3];
    
%Put the barcode/patient id box in the top write corner of the form
% [r c] = size(patient_id_box);
% form(1:r, end-c+1:end) = patient_id_box;
[r c] = size(pma_barcode);
form(1:r, end-c+1:end) = pma_barcode;

%display and save form
figure; imagesc(form); axis image; colormap(gray(256));
imwrite(form, 'C:\isbe\density\form\components\form_image.bmp');