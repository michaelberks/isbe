function [] = rolling_dice_for(n_times)
%Function to simulate rolling a dice n times

%Use a for loop
for i_roll = 1:n_times
    
    %Roll the dice!
    dice_score = roll_a_dice();
    
    %Display the current score
    display(['Roll ' num2str(i_roll) ', you scored ' num2str(dice_score)]);
    
end