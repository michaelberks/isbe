function d = dist(m, n)
%compute the euclidian distance between two N_dimensional vectors
%
% Inputs - m, n : k x N matrices where each row is an N-d vector
%
% Outputs - d: k x 1 vector each row i is the euclidiean distance between
%               the ith rows of m and n

d = sqrt(sum((m-n).^2, 2));