%--------------------------------------------------------------------------
% ------------------------- Week 4 ----------------------------------------
%--------------------------------------------------------------------------

%Script showing various examples of for loops

%--------------------------------------------------------------------------
%%

%1) The most common form of for loop
for iter_var = 1:4
    %Code that does something useful
    display('*** Next iteration in example 1 ***');
    display(iter_var);    
end
display('---------------------');
%%
%2) (1) is equivalent to:
A = 1:4;
for iter_var = A
    %Code that does something useful
    display('*** Next iteration in example 2 ***');
    display(iter_var);
end
display('---------------------');
%%
%3) What if A is a column not a row vector?
A = (1:4)';
for iter_var = A
    %Code that does something useful
    display('*** Next iteration in example 3 ***');
    display(iter_var);    
end
display('---------------------');
%Note 3) has only 1 iteration - because A has only one column, in that
%iteration, iter_var is 4x1 vector
%%
%4) What if A is a matrix?
A = reshape(1:12, 3, 4);
for iter_var = A
    %Code that does something useful
    display('*** Next iteration in example 4 ***');
    display(iter_var);    
end
display('---------------------');
%There are 4 iterations (because A has 3 columns), in each one, iter_var is
%a 3x1 vector
%%
%A doesn't need to contain integers
%5)
A = rand(2,3);
for iter_var = A
    %Code that does something useful
    display('*** Next iteration in example 5 ***');
    display(iter_var);    
end
display('---------------------');

%6)
for iter_var = (0:5)/10;
    %Code that does something useful
    display('*** Next iteration in example 6 ***');
    display(iter_var);    
end
display('---------------------');
%%
% 7) A can be a char array...
A = 'Monday';
for iter_var = A
    %Code that does something useful
    display('*** Next iteration in example 7 ***');
    display(iter_var);    
end
display('---------------------');
%%
%Or a cell array...
%8)
A = {'Monday', 'Tuesday', 'Wednesday'};
for iter_var = A
    %Code that does something useful
    display('*** Next iteration in example 8 ***');
    display(iter_var);    
end
display('---------------------');
%Note each value of iter_var is a 1x1 cell. i.e. iter_var = A(:,i); NOT iter_var = A{:,i};

%9)
A = {'Monday'; 'Tuesday'; 'Wednesday'};
for iter_var = A
    %Code that does something useful
    display('*** Next iteration in example 9 ***');
    display(iter_var);    
end
display('---------------------');
%Note again, the different behaviour because A is now a 3x1 (rather than 1x3) array

%%
%Nested loops (plus a reminder of the flexibility of data types in cells) -
%let's repeat the above examples, using a one big outer loop...
A = cell(9, 1);
A{1} = 1:4;
A{2} = 1:4;
A{3} = (1:4)';
A{4} = reshape(1:12, 3, 4);
A{5} = rand(2,3);
A{6} = (0:5)/10;
A{7} = 'Monday';
A{8} = {'Monday', 'Tuesday', 'Wednesday'};
A{9} = {'Monday'; 'Tuesday'; 'Wednesday'};

for iter_var1 = 1:length(A) %Note the use of length A - this means if I add more examples to A, I don't have to change my code
    %iter_var1 will take values, 1,2,3,...,9 and can index into A
    for iter_var2 = A{iter_var1}
        %Code that does something useful
        display(['*** Next iteration in example ' num2str(iter_var1) ' ***']);
        display(iter_var2);    
    end
    display('---------------------');
end