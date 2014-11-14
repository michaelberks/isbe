function y = nested_example_main_fun1(a, b)
    y = aux_fun(a) + aux_fun(b);
    
end %If end used here 

function y = aux_fun(a)

    y = a^2;
end %All nested functions MUST also have end

% In this example, aux_fun lies outside the main functions end point, thus
% aux_fun has entirely separate memory, thus we can reuse 'a' as a variable

%Aux fun behaves like any other Matlab function - except that only the main
%fun can use it

