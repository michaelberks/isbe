function mam = orientate_mammogram(mam, right)

%Compute the sum of intensities in the 2nd and 3rd quartile of columns
[r c] = size(mam);

r1 = round(r / 3);
r2 = 2*r1;
c_half = round(c / 2);

mam_mean = mean(mam(:));
centre_mean = mean(mam(r1:r2,:));

sum1 = sum(centre_mean(1:c_half) > mam_mean);
sum2 = sum(centre_mean(c_half+1:end) > mam_mean);

%if sum1 > sum2 then the breast is on the left of the image
if (sum1 > sum2) == right
    mam = rot90(mam, 2);
end