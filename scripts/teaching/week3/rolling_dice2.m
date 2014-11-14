function [you_won, roll_count] = rolling_dice2(max_rolls)
%Function to simulate rolling a dice until we score a six - however, in a
%freak of probability, this could go on for ever, and I only want to let
%you win if you roll a 6 within some fixed number of rolls

if nargin < 1
    %Note the example of an if clause, to make the input agrument max_rolls
    %optional, we set a default value here to 3
    max_rolls = 3;
end

%Create a variable to count the number of rolls
roll_count = 0;

%Create a boolean variable, that will be the condition of the loop
keep_rolling = true;

%Create a boolean variable to say if you won
you_won = false;

%Now start a while loop, that will run on the condition dice score isn't 6
while keep_rolling;
    
    %Roll the dice!
    dice_score = ceil(6*rand);
    
    %Increment the roll count
    roll_count = roll_count + 1;
    
    %Display the current score
    display(['Roll ' num2str(roll_count) ', you scored ' num2str(dice_score)]);
    
    %Now check the score and update keep rolling, note the use of if/else
    if dice_score == 6
        %we can stop rolling - so set conditon variable to false
        keep_rolling = 0;
        
        %We won,
        you_won = 1;
        
    elseif roll_count == max_rolls
        %We've used up all our rolls, and not score a six, we have to stop
        %rolling
        keep_rolling = 0;
    end %Otherwise we just keep going, so don't need an else
    
end

%Did we win?
if you_won
    display(['Congratulations, you won after ' num2str(roll_count) ' rolls!!']);
else
    display(['I''m sorry, you used up ' num2str(max_rolls) ' rolls without scoring a 6. You lose!!']);
end