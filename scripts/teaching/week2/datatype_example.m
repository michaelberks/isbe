%Script showing examples of different data types. Also note the use of cell
%mode so you can execute batches of commands together
%% Numerical data types
% Doubles, uint8's

%By default, arrays are created with double type
x_double = zeros(4);

%But can 'cast' to other types
x_uint8 = uint8(x_double);

%Use the function whos at the command line, to see details of the variables
%in memory
whos; %Note how x_uint8 use 1/8 of the memory of x_double

%Can also use ones and doubles to create data types other than double
%directly (note the use of an optional argument - we'll come back to this later!)
y_uint8 = zeros(4, 'uint8');
z_uint8 = ones(4, 'uint8');
whos;
clear;
%% Combining different types
% Different types in general can't be defined
x_double + x_uint8; %generates error - trying to combine a double and an integer type

%Some operations and functions aren't defined for some types
x_uint8^2; %generates error - matrix multiplication is not defined for integer classes
x_uint8(1)^2; %ok - because, a single element uses scalar multiplication

%Bounds are automatically capped
z_low = uint8(-1);
z_high = uint8(256);

%Decimal values are automatically rounded to the nearest integer (in range)
r_double = rand(4);
r_uint8 = uint8(r_double);
clear;
%% Characters and strings
a = 'a'; %A single char
abc = 'abcdefghi'; %A row vector of chars forms a string
my_string = 'Happy Birthday!!'; %And another

display(abc(5)); %Can index like a numeric vector
display(my_string(1:5));

abc2 = reshape(abc, 3, 3); %Can also create 2D arrays

%Can create strings from numeric types using num2str - useful for display
%output, creating filenames etc.
r = rand(1);
r_char = num2str(r);
size(r);
size(r_char);

%Can combine char arrays just like numeric arrays
for i_im = 1:5 %More on for loops next week!
    im_name = ['my_image_' num2str(i_im) '.jpg'];
    display(im_name);
    %im = imread(im_name); for example...
end
clear;
%% Logical data - true or false?
X = rand(5, 3); display(X);
b = X < 0.5; display(b); %b is 5x3 array of logical values, 1 (true) where X is less than 0.5, 0 (false) otherwise
y = X(b); display(y); %y is a column vector of all the values in X less than 0.5

%Can create directly
i_double = zeros(5, 3);
i_f = false(5, 3);
i_t = true(5, 3);

%Can index the logical arrays just like doubles
i_double(2,2) = 1;
i_f(2,2) = 1;
i_t(2,2) = 0;
whos;

%Can use logical arrays to index any other arrays - returns only those
%parts of the array where the logical index is 1
y_f = X(i_f); display(y_f); %The 2nd row of the 2nd column in X
y_t = X(i_t); display(y_t); %Everything except the 2nd row of the 2nd column in X
y_d = X(i_double); %Generates error - although i_double contains only 0s and 1s, it must be stored as a logical to index in this way

