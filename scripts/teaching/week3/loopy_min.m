function [min_val, min_row, min_col] = loopy_min(array_in)
%What if we didn't have the min function and wanted to work out the minimum value if an array?

%Get size of array
[rows cols] = size(array_in);

%Set up a variable to store the min value
min_val = inf;

%Now loop through the array along the rows and columns
for i_row = 1:rows
    for i_col = 1:cols
        
        %check each value in input array against the current min_val
        if array_in(i_row,i_col) < min_val
            
            %if it's smaller, update the min val and row, col positions
            min_val = array_in(i_row,i_col);
            min_row = i_row;
            min_col = i_col;
        end
    end
end

%If array_in has several copies of the same minimum value, which one will
%this return? How might we change this?

%What happens if array_in isn't a numerical 2D array? How can we check for
%this?


