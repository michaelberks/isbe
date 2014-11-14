function y = my_fun2(a, b)
%My second function, uses MY_FUN to generate a new quadratic function in a
%and b
 y = my_fun(a,b) + my_fun(2*a, b/2);
