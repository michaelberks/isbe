function focus_expt()

clc;

% imgroot = 'U:\projects\nailfold\capture\2013_03_06\Left\Digit4\x300\';
% imgpath = fullfile(imgroot,'autofocus_pattern');
% imgpath = fullfile(imgroot,'autofocus_vessels');

imgroot = 'U:\projects\nailfold\capture\2013_03_07\Left\Digit4\x300\';
imgpath = fullfile(imgroot,'17_56_56');

% dat = dlmread(fullfile(imgpath,'sharpness.log'));
% dat = dat(:,end-1:end);
% figure(2); clf; hold on;
%     plot(dat(:,1), dat(:,2), 'b.');
%     ylim([0,30]);

d = dir(fullfile(imgpath, 'frame_*.png'));
f_inds = 1:1:min(251, length(d));
n = length(f_inds);

d = d(f_inds);
sharpness = nan(n, 8);
leg = repmat({''}, [1, size(sharpness, 2)]);

for i = 1:length(d)
    img = double(imread(fullfile(imgpath,d(i).name)));
    img = img - mean(img(:));
%     img = img / (std(img(:)) / 40);
    img = img + 128;
%     img = normim(img, 'stretch') * 255;
    
    sharpness(i,1) = sharpness_sobel(img); leg{1} = 'sobel';
    sharpness(i,2) = img_var(img); leg{2} = 'var';
    
    roi = get_roi(img);
    sharpness(i,3) = var(roi(:)); leg{3} = 'totvar';
%     sharpness(i,3) = img_skew(img); leg{3} = 'skew';
%     sharpness(i,2) = img_entropy(img); leg{2} = 'entropy';
%     sharpness(i,3) = img_range(img); leg{3} = 'range';
%     sharpness(i,5) = img_dog(img); leg{5} = 'DoG';
%     sharpness(i,6) = img_g1d(img); leg{6} = 'g1d';
%     sharpness(i,7) = img_corr(img); leg{7} = 'corr';

%     [xx,yy] = meshgrid(1:size(img,2), 1:size(img,1));
%     figure(10);
%         plot3(xx(:), yy(:), img(:), 'b.', 'markersize', 1);

    disp([f_inds(i) sharpness(i,:)]);
    
%     pause;
end
% 
% sharpness = sharpness - repmat(min(sharpness), [n, 1]);
sharpness = sharpness ./ repmat(max(sharpness), [n, 1]);

figure(3); clf; hold on;
    plot(f_inds, sharpness, '.-');
    legend(leg);
ylim([-0.1, 1.1]);
% ylim([0,1]);
    

function sharpness = sharpness_sobel(img)

sobel_x = conv2([1 2 1]/4, [1 0 -1]/2, img, 'same');
sobel_y = conv2([1 0 -1]/2, [1 2 1]/4, img, 'same');
sobel_mag = sqrt(sobel_x.^2 + sobel_y.^2);

sobel_mag = get_roi(sobel_mag);

map = (sobel_mag > 0);
sharpness = mean(sobel_mag(map));


function value = img_corr(img)

img00 = img(1:2:end, 1:2:end);
img10 = img(1:2:end, 2:2:end);
img01 = img(2:2:end, 1:2:end);
% img11 = img(2:2:end, 2:2:end);

a = [img00(:); img00(:)]; 
a = a - mean(a);

b = [img01(:); img10(:)]; 
b = b - mean(b);

figure(5); clf; hold on;
    plot(a(:), b(:), 'b.', 'markersize', 1);
    axis('equal', 50*[-1,1,-1,1]);
    
value = (a'*b) / sqrt((a'*a)*(b'*b));


function value = img_g1d(img)

[g,dg,ddg] = gaussian_filters_1d(5);
dg_x = conv2(g, dg, img, 'same');
dg_y = conv2(dg, g, img, 'same');
dg_mag = sqrt(dg_x.^2 + dg_y.^2);

roi = get_roi(dg_mag);

map = (roi > 0);
value = mean(roi(map));


function value = img_range(img)

roi = get_roi(img);
value = max(roi(:)) - min(roi(:));
value = double(value);


function value = img_entropy(img)
roi = get_roi(img);
h = hist(roi(:), 0:255);

p = h / sum(h(:));

value = -p(p>0) * log(p(p>0))';

% figure(1); clf;
%     bar(p);
%     xlim([120, 160]); ylim([0, 0.15]);


function value = img_skew(img)
roi = get_roi(img);

mn = mean(roi(:));
sd = std(roi(:));

value = mean( ((roi(:)-mn) ./ sd(:)).^3 );
value = abs(value);

% figure(1); clf; colormap(gray(256));
%     imagesc(roi);


function value = img_dog(img)
dog = [ 1  2   1;
        2 -12  2;
        1  2   1] / 16;
img = conv2(img, dog, 'same');
roi = get_roi(img);

% figure(1); clf; colormap(gray(256));
%     imagesc(roi);
%     pause;
    
value = mean(roi(:));


function value = img_var(img)
img = get_roi(img);
roi = [zeros(size(img,1), 1) img];
roi = [zeros(1, size(roi,2)); roi];

intimg = cumsum(cumsum(roi, 1), 2);
intimg2 = cumsum(cumsum(roi.^2, 1), 2);

hw = 10;
fw = 2*hw+1;
fw2 = fw*fw;

imgvar = nan(size(roi));
for i = 2+hw:size(roi, 2)-hw
    for j = 2+hw:size(roi, 1)-hw
        sm = intimg(j+hw,   i+hw) + intimg(j-hw-1, i-hw-1) - ...
             intimg(j-hw-1, i+hw) - intimg(j+hw,   i-hw-1);
         
        sm2 = intimg2(j+hw,   i+hw) + intimg2(j-hw-1, i-hw-1) - ...
              intimg2(j-hw-1, i+hw) - intimg2(j+hw,   i-hw-1);
         
        imgvar(j, i) = (sm2 / fw2) - (sm / fw2)^2;
    end
end

imgvar = imgvar(2+hw:end-hw, 2+hw:end-hw);
map = (imgvar > 0);
% value = mean(imgvar(map)) / var(img(:));
value = mean(imgvar(map));

imgvar(~map) = 0;
% figure(10); clf; colormap(gray(256));
% 	imagesc(imgvar);
% figure(5); clf;
%     hist(imgvar(:), 0:50);
    
    
function roi = get_roi(img)
ni = round(size(img, 2) / 2);
nj = round(size(img, 1) / 2);

i0 = round((size(img, 2) - ni) / 2);
j0 = round((size(img, 1) - nj) / 2);

step = 1;
for i = 1:log2(step)
    img = conv2([1,2,1]/4, [1,2,1]/4, img, 'same');
end

roi = img(j0:step:j0+nj-1, i0:step:i0+ni-1);
