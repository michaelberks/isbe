function [step_g, errorcheck] = calculate_step_values(mammo, step_centres_xy, debug_mode)
%CALCULATE_STEP_VALUES Compute mean grey level in the centre of each step
%   [step_g, errorcheck] = calculate_step_values(mammo,step_centres_xy,debug_mode)
%
% Inputs:
%      mammo - Mammogram containing stepwedge
%
%      step_centres_xy - Coordinates of each step centres
%
%      debug_mode - flag to determine debug output
%
%
% Outputs:
%      step_g - Mena grey-level at each step centre
%
%      errorcheck - Flag if anything went wrong
%
%
% Example:
%
% Notes: Only require user input in debug mode (even in the case of errors)
%
% See also:
%
% Created: 25-Sep-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester
if nargin < 3
    debug_mode = 0;
end

errorcheck = 0;
num_steps = size(step_centres_xy,1);
step_g = zeros(num_steps,1);

%Set default step sizes (i.e. the area of the step)
step_l = 3;
step_w = 5;

for ii = 1:num_steps

    % define edges of this step's region
    xmin = round(step_centres_xy(ii,1)-step_w);
    xmax = round(step_centres_xy(ii,1)+step_w);
    ymin = round(step_centres_xy(ii,2)-step_l);
    ymax = round(step_centres_xy(ii,2)+step_l);

    %extract step from mammogram
    this_step = mammo(ymin:ymax,xmin:xmax);

    step_g(ii) = mean(double(this_step(:)));

    %Check standard deviation of step, probably large if region doesn't lie
    %fully within one step
    if debug_mode
        wedge_std = std(double(this_step(:)));
        if(wedge_std > 1000)
            if debug_mode
                disp('step has large std - look at intensity histogram');
                disp('press enter to continue...');
                pause;
                figure;
                imhist(this_step);
            end
        end
    end
end

%Check we've got real values for every step
if any(isnan(step_g)) 
    if debug_mode
        repeat_step = 'Repeat stepwedge selection';
        skip = 'Skip this mammogram';
        answer = questdlg('An error has occurred selecting step values. Please either:',...
            'Stepwedge error', repeat_step, skip, repeat_step);

        if strcmpi(answer, skip);
            errorcheck = 1;
            return;
        else
            [step_centres_xy] = locate_stepwedge_manual(mammo, debug_mode);
            [step_g, errorcheck] = calculate_step_values(mammo, step_centres_xy, debug_mode);
            return
        end
    else
        errorcheck = 1;
        return;
    end
end

%--------------------------------------------------------------------------
% Now make sure the steps have monotonically increasing values. For steps
% 1-20, correct the lower step. For steps 21-39, correct the upper step

%If we're debugging, plot the step values before/after any corrections
if debug_mode
    figure('Name', 'step values'); hold on;
    plot(step_g,'bo');
    xlabel('step number');
    ylabel('mean intensity');
    title('stepwedge intensity vs step number');
end

halfwaystep = round((num_steps-1)/2);
for ii = halfwaystep:num_steps-1
    if(step_g(ii) < step_g(ii+1))
        step_g(ii+1) = step_g(ii);
        
        if debug_mode
            disp(['step above ', num2str(ii), ' has intensity higher - setting same']);
        end
        
    end
end
for ii = halfwaystep:-1:2
    if(step_g(ii-1) < step_g(ii))
        step_g(ii-1) = step_g(ii);
        
        if debug_mode
            disp(['step below ', num2str(ii), ' has intensity higher - setting same']);
        end
        
    end
end

%Plot the new values following any corrections
if debug_mode
    plot(step_g,'r+');
    legend('before correction','after correction');
end
