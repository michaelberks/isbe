%--------------------------------------------------------------------------
% ------------------------- Week 4 ----------------------------------------
%--------------------------------------------------------------------------

%Script showing various examples of if/else statments

%--------------------------------------------------------------------------
%% 1) A simple if statement to display if a number is even

a = ceil(rand*10); %Generates a random integer between 1 and 10
if rem(a, 2) == 0
    display([num2str(a) ' is even.']);
end

%% 2) A simple if statement to display if a number is even or odd

a = ceil(rand*10); %Generates a random integer between 1 and 10
if rem(a, 2) == 0
    display([num2str(a) ' is even.']);
else
    display([num2str(a) ' is odd.']);
end

%% 3) Using elseif 
a = ceil(rand*10) / 2;
if rem(a, 2) == 0
    display([num2str(a) ' is even.']);
elseif rem(a, 2) == 1
    display([num2str(a) ' is odd.']);
else
    display([num2str(a) ' is not an integer.']);
end

%% 4) Using nested clauses, instead of elseif
% (4) is identical to (3)
if rem(a, 2) == 0
    display([num2str(a) ' is even.']);
else
    if rem(a, 2) == 1
        display([num2str(a) ' is odd.']);
    else
        display([num2str(a) ' is not an integer.']);
    end
end

%% 5) Using a switch statement instead of if/else
%This does exactly the same as 3 and 4
a = ceil(rand*10) / 2;
rem_val = rem(a, 2);
switch rem_val
    case 0
        display([num2str(a) ' is even.']);
        
    case 1
        display([num2str(a) ' is odd.']);
        
    otherwise
        display([num2str(a) ' is not an integer.']);
end