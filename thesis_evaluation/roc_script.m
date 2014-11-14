% ROC script

% Some test code for computing ROC curves for each observer in an observer
% trial. To be expanded into to proper functions etc.

% Each observer rates each mass on a 5 point scale:
% 1) Definitely synthetic
% 2) Probably synthetic
% 3) Possibly real
% 4) Probably real
% 5) Definitely real

% Let's run some repeat trials
figure('Name','ROC Curve for random assignment of ratings'); 
for rp = 1:12
    
    % Generate some random ratings for the real masses:
    real_ratings = randsample(5, 30, true);

    %Generate some random weighting for the synthetic masses
    synthetic_ratings = randsample(5, 30, true);

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

    subplot(3,4,rp); plot(roc_curve(:,1), roc_curve(:,2)); axis equal; axis([0 1 0 1]);
    AUC = sum( (roc_curve(1:5,1)-roc_curve(2:6,1)) .* roc_curve(2:6,2)) + ...
        0.5*sum( (roc_curve(1:5,1)-roc_curve(2:6,1)) .* (roc_curve(1:5,2)-roc_curve(2:6,2)) );
    text(0.5, 0.4, ['AUC = ', num2str(AUC)]);
end

%%
%Let's do that again but now weight the sampling to show differentiation
weights_real = (1:5) / 15;
weights_synthetic = (5:-1:1) / 15;
figure('Name','ROC Curve for random assignment of ratings'); 
for rp = 1:12
    
    % Generate some random ratings for the real masses:
    real_ratings = randsample(5, 30, true, weights_real);

    %Generate some random weighting for the synthetic masses
    synthetic_ratings = randsample(5, 30, true, weights_synthetic);

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

    subplot(3,4,rp); plot(roc_curve(:,1), roc_curve(:,2)); axis equal; axis([0 1 0 1]);
    AUC = sum( (roc_curve(1:5,1)-roc_curve(2:6,1)) .* roc_curve(2:6,2)) + ...
        0.5*sum( (roc_curve(1:5,1)-roc_curve(2:6,1)) .* (roc_curve(1:5,2)-roc_curve(2:6,2)) );
    text(0.5, 0.4, ['AUC = ', num2str(AUC)]);
end
%%
aucs = zeros(10,1);
for ii = 1:10
    load(['C:\isbe\dev\observer_study\user_data\user', zerostr(ii,2)]);
    
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
    plot(roc_curve(:,1), roc_curve(:,2)); axis equal; axis([0 1 0 1]); hold on;
    plot(roc_curve(1,1), roc_curve(1,2), 'rx'); axis equal;
    AUC = sum( (roc_curve(1:5,1)-roc_curve(2:6,1)) .* roc_curve(2:6,2)) + ...
        0.5*sum( (roc_curve(1:5,1)-roc_curve(2:6,1)) .* (roc_curve(1:5,2)-roc_curve(2:6,2)) );
    text(0.5, 0.4, ['AUC = ', num2str(AUC)]);
    aucs(ii) = AUC;
    clear user_data;
end
%%
%Do ROC again using calculate_roc_curve function
aucs = zeros(10,1);
times = zeros(10,1);
for ii = 1:10
    load(['C:\isbe\dev\observer_study\user_data\user', zerostr(ii,2)]);
    
    class_labels = user_data.ratings(:,1) <= 30;
    class_values = user_data.ratings(:,2);
    [roc_curve auc] = calculate_roc_curve(class_values,class_labels,0:5);
    

    figure;
    plot(roc_curve(:,1), roc_curve(:,2)); axis equal; axis([0 1 0 1]); hold on;
    plot(roc_curve(1,1), roc_curve(1,2), 'rx'); axis equal;

    text(0.5, 0.4, ['AUC = ', num2str(auc)]);
    aucs(ii) = auc;
    times(ii) = mean(user_data.ratings(:,3));
    clear user_data;
end