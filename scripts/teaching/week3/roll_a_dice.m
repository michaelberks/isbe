function [score] = roll_a_dice()
%Function to simulate rolling a dice - i.e. returns a random integer
%between 1 and 6

%Roll the dice!
score = ceil(6*rand);