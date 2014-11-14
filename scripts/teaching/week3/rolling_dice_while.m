function [roll_count] = rolling_dice_while()
%Function to simulate rolling a dice until we score a six

%Create variable to store the current dice score, set to zero before we've
%thrown for the first time
dice_score = 0;

%Create a variable to count the number of rolls
roll_count = 0;

%Now start a while loop, that will run on the condition dice score isn't 6
while dice_score ~= 6;
    
    %Roll the dice!
    dice_score = roll_a_dice();
    
    %Increment the roll count
    roll_count = roll_count + 1;
    
    %Display the current score
    display(['Roll ' num2str(roll_count) ', you scored ' num2str(dice_score)]);
    
end

%The loop above continue until dice score equals 6, then control returns to
%the main program and the following lines are executed
display(['Congratulations, you won after ' num2str(roll_count) ' rolls!!']);