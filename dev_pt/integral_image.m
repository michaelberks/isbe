function intimg = integral_image(img)
% function to compute the integral image (summed area table)

intimg = cumsum(cumsum(img,1),2);

