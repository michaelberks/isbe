function [] = roc_gui()
%ROC_GUI *Insert a one line summary here*
%   [] = roc_gui()
%
% Inputs:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 06-Jul-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

go_on = true;
user_path = 'C:\isbe\dev\observer_study\user_data\';

while go_on
    
    [user_file, user_path] =...
        uigetfile('*.mat','Select the user data to process', user_path, 'Multiselect','off');
    
    if ~user_file
        answer = questdlg(...
            'Select another file?','No file selected','Yes', 'Exit', 'Yes');
        if strcmpi(answer, 'yes')
            user_path = 'C:\isbe\dev\observer_study\user_data\';
            continue;
        else
            go_on = false;
            continue;
        end
    end
    
    temp = load([user_path, user_file]);
    user_data = temp.user_data; clear temp;
    
    real_mass_idx = user_data.ratings(:,1) <= 30;
    syn_mass_idx = user_data.ratings(:,1) > 30;
    real_ratings = user_data.ratings(real_mass_idx,2);
    synthetic_ratings = user_data.ratings(syn_mass_idx,2);

    roc_curve = zeros(6,2);

    % Calcultae ROC points at thresholds from 0 to 5
    for th = 0:5
        

        %Compute true positives, and true positive fraction
        TP = sum(real_ratings > th); %also FN = 30 - TP;
        TPF = TP / 30; %since TPF = TP/(TP+FN) = TP/30;

        %Compute false positives, and false positive fraction
        FP = sum(synthetic_ratings > th); %also TN = 30 - FP;
        FPF = FP / 30; %since FPF = FP/(FP+TN) = TP/30;

        roc_curve(th+1,1) = FPF;
        roc_curve(th+1,2) = TPF;

    end

    figure;
    plot(roc_curve(:,1), roc_curve(:,2), 'LineWidth', 2); 
    axis equal; axis([-0.01 1.01 -0.01 1.01]);
    AUC = sum( (roc_curve(1:5,1)-roc_curve(2:6,1)) .* roc_curve(2:6,2)) + ...
        0.5*sum( (roc_curve(1:5,1)-roc_curve(2:6,1)) .* (roc_curve(1:5,2)-roc_curve(2:6,2)) );
    text(0.5, 0.4, ['AUC = ', num2str(AUC)]);
    title(['ROC analysis curve for user ', num2str(str2num(user_file(5:6)))]); %#ok
    xlabel('False positive fraction');
    ylabel('True positive fraction');
    clear user_data;
    
    answer = questdlg(...
        'Select another file?','ROC analysis complete','Yes', 'Exit', 'Yes');
    
    if strcmpi(answer, 'exit')
        go_on = false;
    end
end
    
    
