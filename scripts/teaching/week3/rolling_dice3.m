function [you_won, roll_count] = rolling_dice3(max_rolls)
%Function to simulate rolling a dice until we score a six - however, in a
%freak of probability, this could go on for ever, and I only want to let
%you win if you roll a 6 within some fixed number of rolls. This works the
%same as rolling_dice2, but uses a for loop and a break statement instead

if nargin < 1
    %Note the example of an if clause, to make the input agrument max_rolls
    %optional, we set a default value here to 3
    max_rolls = 3;
end

%Create a boolean variable, that will be the condition of the loop
keep_rolling = true;

%Create a boolean variable to say if you won
you_won = false;

%This for loop will run from 1 to max_rolls (unless we trigger the break
%statement)
for roll_count = 1:max_rolls
    
    %Roll the dice!
    dice_score = ceil(6*rand);
    
    %Display the current score
    display(['Roll ' num2str(roll_count) ', you scored ' num2str(dice_score)]);
    
    %Now check the score
    if dice_score == 6
        %We won, and can stop rolling by breaking out the loop
        you_won = 1;
        break;
    end    
end

%Did we win?
if you_won
    display(['Congratulations, you won after ' num2str(roll_count) ' rolls!!']);
else
    display(['I''m sorry, you used up ' num2str(max_rolls) ' rolls without scoring a 6. You lose!!']);
end