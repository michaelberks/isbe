clc;

N = 4;
X = rand(N,3);
beta = rand(3,1);
y = X*beta;

X2 = [X ones(N,1)];
y2 = y+1;
beta2 = X2\y2

X2*beta2 - y2