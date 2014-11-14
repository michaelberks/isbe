function y = nested_example_main_fun3(a, b)
    y = aux_fun(2) + b;
    
 %We can also omit the end, but... 

function y = aux_fun(c)

    y = a^c;
end
%The nested function now lies inside the main function. This means aux_fun
%shares the main function's memory. Any variable in the main function can
%also be called (and modified) in the nested function. In this case, we
%have shared a - in more recent versions of Matlab the editor depicts 'a'
%with a different colour to tell us we've shared it across multiple
%functions

%Note that this can get confusing, and should generally be avoided -
%however there are some occasions (you may come across one in Tim's
%optimisation exercises) where it is a useful feature to have


end %Using nesting in this form, we MUST include an end for the main function


