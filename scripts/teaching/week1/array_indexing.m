%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Example Script: Indexing arrays
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Using indexing to access parts of an arrray
%--------------------------------------------------------------------------
A = reshape(1:16, 4, 4); %Create a 4x4 matrix with elements 1 to 16
display(A);

%Access individual elements by giving row and column subscripts or giving a single index
b = A(2, 3) %3rd column  of the 2nd row
b = A(10) %10th element - also the 3rd column of the second row

%Access whole rows or columns using :
b = A(2,:) %b is 1x4, the 2nd row
c = A(:,4) %c is 4x1, the 4th column

%Access specific ranges
b = A(1, 2:4) %b is 1x3, the 2nd, 3rd and 4th elements of the 1st row
C = A(2:3,2:3) %C is 2x2, the middle 4 elements of A
D = A(2:3, [1 4]) %D is 2x2, the 1st and 4th elements of the 2nd and 3rd rows

%Using the keyword 'end'
A = rand(4,5)
b = A(1, 2:end) %b is 1x4, the 2nd to 5th elements of the 1st row
c = A(2:end-1,3) %c is 2x1, the 3rd column of the 2nd to 3rd rows
D = A(1:end/2,:) %D is 2x5, all columns of the first 2 rows

%--------------------------------------------------------------------------
%Using indexing to assign values to an arrray
%--------------------------------------------------------------------------
%Create A and B
A = zeros(5,5) 
B = rand(4,4)

%Assign parts of A using parts of B
A(1:4,1:4) = B %Copies all of B to the top-left of A
A(3:5,4:5) = B(1:3,1:2) %Copies the top-left of B to the bottom-right of A
A(3:5,4:end) = B(1:3, [1 4]) %Copies the first and last columns of the first 3 rows of B to to the bottom right of A
A(1:4,1:3) = B %Causes an error! LHS is 4x3, RHS is 4x4

A([1 5 9]) = B([2 4 11]) %Copies the 2nd, 4th and 11th elements of B to the 1st, 5th and 9th elements of A
A([3 7 8]) = B([2 11 2]) %Copies the 2nd and 11th elements of B to the 3rd, 7th and 8th elements of A - note B(2) is copied to both A(3) and A(8)
A([3 7 8]) = B([2 11]) %Causes an error! LHS is 1x3, RHS is 1x2

%Copy back to the same variable
B(1,:) = B(4,:) %Copies the 4th row of B to its 1st row (the values of which are now lost!)

%Assign with a scalar
A(:,2) = 0 %Sets the 2nd column of A to zero
B([1 4 15]) = pi %Sets the 1st, 4th and 15th elements of B to ?

%-------------------------------------------------------------------------
%Out of bounds indexing - access
%-------------------------------------------------------------------------
A = rand(5,5)
b = A(5,6) %Trying to access the 6th column
b = A(26) %Trying to access the 26th element
b = A(end+1,2) %Trying to access to 6th row

%-------------------------------------------------------------------------
%Out of bounds indexing - assigning
%-------------------------------------------------------------------------
A(5,6) = 1 %A is now 5x6 - note how the remaining elements in the 6th column are automatically set to zero
A(end+2,2) = 2 %A is now 7x6
A(:,end+1) = rand(7,1) %A is now 7x7
A(-1,1) = 0 %Causes error, indices can't be negative
A(1.5,1) %Causes error, indices must be integer values
A(1:2.5,1) %Causes a warning not an error, why?