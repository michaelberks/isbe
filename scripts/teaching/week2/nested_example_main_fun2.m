function y = nested_example_main_fun2(a, b)
    y = aux_fun(a) + aux_fun(b);
    
%end %We can also omit the end, but... 

function y = aux_fun(a)

    y = a^2;
%end %We must also then omit the end from the nested functions

%aux fun behaves in the same way as in example, i.e. it is an entirely
%separate function


