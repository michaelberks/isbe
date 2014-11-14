%**************************************************************************
%A matlab tutorial script - some basic functions messing around with images
%as well other simple introductions to matlab
%**************************************************************************

%Comments: comments are created using the '%' character, any text after the
% '%' will be interpreted as a comment

%Help: documentation for any of the core functions used below can be
%obtained by calling >>help fun_name at the command line or through the
%interactive help menu

%Cell mode, this script uses cell mode which may not be available on the
%msc machines (dpending on the version of matlab they have installed).
%Cells are demarcated by '%%' characters and can be executed individually
%using the buttons on the toolbar

%% - loading in an image
% load in one of matlab's stock images using the imread function - the image
% will be stored in an array with data type depending on the image format.
% Color images will be loaded as unsigned 8 bit integer arrays of size [r x c x 3] 
% arrays coding the RGB channels.
% Grayscale images will be loaded as unsigned 8 or 16 bit [r x c] arrays
im = imread('cameraman.tif'); 

%Get dimensions of image
[r c] = size(im); %note that the ';' suppresses output
[r2 c2] = size(im)  %without the ';' the result of the command will be displayed in the command window

%Display the image - for now we'll just use imshow although there's more
%discussion of images in section ?
imshow(im);

%% - basic processing of arrays
%Now an image is stored as array, we can apply all sorts of operations to
%it using the core matlab functions

%Compute the mean intensity of the image
m = mean(im(:));
%Note we use im(:) - this collapses the array into a single column vector -
%like nearly all array operation in matlab, mean by default operates on the
%columns of 2D array. So for example, the command below computes the mean
%of each column
m_col = mean(im); %m_col will be a 1xc row vector;

%we can specify a dimension as a second argument to mean, so to compute the
%mean of each row
m_row = mean(im, 2); %m_row will be rx1 column vector

% Elements of an array are accessed using '(' brackets in r, c format. So
% im(i, j) accesses the element in the j-th column of i-th row.
% Alternatively, each element can be accessed by a single index that
% counts column wise, so if k = r*(j-1) + i, then im(k) = im(i,j)
%
% ':' on its own specifies that you select all the elements along that dimension, so
% for example im(i,:) references the i-th row of im

%So to access the first row of im
r1 = im(1,:);

%However, when placed between start and finish values, ':' generates a
%vector of values e.g.
1:10 %creates the row vector [1 2 3 4 5 6 7 8 9 10]
(1:10)'; %Uses the transpose operator to form a column vector

%You can also specify a different step size other than 1, e.g.
1:2:10 %creates the row vector [1 3 5 7 9]

%This can be used to index arrays, so to obtain the first 10 elements in
%the first row of im
r1_1to10 = im(1,1:10);

%Or the last 10
r1_last10 = im(1,c-9:c);

%Or every odd/even element in the first row
r1_odd = im(1,1:2:c);
r1_even = im(1,2:2:c);

%When used within the indexing brackets of any array 'end' acts as a
%keyword for the last element, so the above command could also be generated
%by... 
r1_last10 = im(1,end-9:end);

%Arrays can be indexed in both dimensions, so to obtain a block of elements
block = im(5:10, 2:3); %returns a 6x2 array

%Of course, providing it contains integer elements, any vector can be used
%to index an array so e.g.
idx = [1:3 7 19]; im(idx,3); %references the 1st, 2nd, 3rd, 7th and 19th elements of the 3 column

%We can alos flip things round and use the pixel values in im, as indices.
% That the elements in im are of the range 0-255, so imagine we have a
% lookup table that assigns some value to each number from 0-255 (e.g. see
% colormaps in displaying images topic below)
lookup = rand(255,1); %creates 256x1 vector of random numbers ~U(0,1)
converted_image = lookup(im); %produces an r x c array of values from the lookup table

% The above operation is extactly the kind of operation we use when
% converting grey level pixel intensities in the digitised mammograms into
% glandular thicknesses at each location in the breat (the stepwedge
% program is essentially about how we compute the lookup table from the
% calibration data in the image)

%% Image masks
%We can also index arrays with logical vectors/arrays e.g.
mask1 = im > 100; %generates a binary r x c map, with 1's at all pixels with intensity greater than 100
N = sum(mask1(:)); %Count how many elements are greater than 100
greater_than100 = im(mask1); % generates an Nx1 vector, of all the pixel intensities greater than 100
less_than100 = im(~mask1); %uses the inverse of the mask to select elements less than or equal to 100

%We can do this without explicity creating the mask as a variable
greater_than50 = im(im > 50);

%We can also operate individually along columns or rows
greater_than75 = im(im(:,1) > 75,:); %references the elements of each column in im for which the values in the first column are greater than 75

%It may also be useful to find out where these pixels are in the image,
%which we can do using find
[mask_idx] = find(mask1); %this produces indices, alternatively
[mask_rows mask_cols] = find(mask1); %produces row/col subscripts

%Note we can always swap between subscripts and indices using...
mask_idx = sub2ind([r c], mask_rows, mask_cols);
[mask_rows mask_cols] = ind2sub([r c], mask_idx); %Note the example of an overloaded function

%We can also create a mask of a geometrical area in an image - in mammograms
% for example, we often use a mask specifying pixels that actually lie in
% the breast
%Below we create a circular mask and use it in conjuntion with the lookup
%operation from above
xx = repmat(1:c, r, 1); %use repmat to generate arrays of x/y coordinates
yy = repmat((1:r)', 1, c);
mask2 = (xx - c/2).^2 + (yy - r/2).^2 < 128^2; %See combining arrays below

converted_image2 = lookup(im);
converted_image2(~mask2) = 0; %Discard all values from the converted image outisde of the mask
figure; imshow(converted_image);
%% Combining arrays
%As in creating circular masks above, +, *, ^ etc can be used to perform
%element-wise arithmetic
im2 = im + im; %elementwise addition - size of arrays must match
im_plus_10 = im + 10; %adding a scalar to an array
im_squared = im .* im; %element wise multiplication - size of arrays must match

%To perform matrix multiplication, omit the '.' - although this isn't
%supported for the uint8 data type (see note on data types below) 
mat_squared = im * im; %will produce error
mat_squared = double(im) * double(im); %matrix multiplication

%You can't automatically apply row/cols to an array though. (e.g. as if
%subtracting the mean from sample data), so
im - im(1,:); %try subtracting the first row from all rows - this will error
im - repmat(im(1,:), r, 1); % will work 
bsxfun(@minus, im, im(1,:));% or use bsxfun (I can never remember its blooody name - stands for binary singleton expansion I think??)

%A quick note on data types: most matlab functions aren't fussy about
%datatypes (uint8, uint16, double etc) but datatypes can't usually be
%combined, so
im_double = double(im);
im + im_double; %will produce error

%Also binary masks need to be logical, otherwise the zeros are interpreted
%as an out-of-bounds index. e.g.
mask3 = zeros(r, c);
im(mask3); %will error
im(logical(mask3)); %is ok

%**************************************************************************
% That's all for now! When I get time I'll add in the topics below
%**************************************************************************

%% - displaying an image: taking more manual control than imshow
%% - reading and writing data
%% - structures, cells
%% - working with multidimensional arrays
%% - image filtering and convolutions
%% - using optional arguments in your own functions





