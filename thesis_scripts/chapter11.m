%
%Figures for the evaluation chapter
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
mkdir C:\isbe\thesis\figures\11
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Compute average score for each mass in the test
num_users = 10;
all_ratings = NaN(90, num_users);
all_timings = NaN(90, num_users);
all_ids = cell(1,num_users);
all_feedback = zeros(90, num_users, 7);
for ii = 1:num_users
    user_id = ['User' zerostr(ii,2)];
    load(['C:\isbe\dev\observer_study\user_data\user', zerostr(ii,2)]);
    for jj = 1:60;
        all_ratings(user_data.ratings(jj,1),ii) = user_data.ratings(jj,2);
        all_timings(user_data.ratings(jj,1),ii) = user_data.ratings(jj,3);
    end
    for jj = 1:length(user_data.feedback_other)
        for kk = 1:7
            all_feedback(user_data.feedback(jj,1),ii,kk) =...
                all_feedback(user_data.feedback(jj,1),ii,kk) + user_data.feedback(jj,kk+1);
        end
    end
    all_ids{ii} = user_id;    
end
avg_score = nanmean(all_ratings,2);
%
%Display averag scores and significance test
display([' Mean of real masses = ', num2str(mean(avg_score(1:30)))]);
display([' Mean of synthetic masses = ', num2str(mean(avg_score(31:end)))]);
display([' Mean % real masses = ', num2str(mean(avg_score(1:30))/5)]);
display([' Mean % synthetic masses = ', num2str(mean(avg_score(31:end))/5)]);
[p,h,stats] = ranksum(avg_score(1:30), avg_score(31:end));
display(['MAIN TEST: Mann-whitney U test, p = ', num2str(p)]);
display(stats);
%%
%--------------------------------------------------------------------------
% Get reader ratings for observer experiment and draw ROC curves
aucs = zeros(10,1);
times = zeros(10,1);

for ii = 1:10
    load(['C:\isbe\dev\observer_study\user_data\user', zerostr(ii,2)]);
    
    class_labels = user_data.ratings(:,1) <= 30;
    class_values = user_data.ratings(:,2);
    [roc_curve auc] = calculate_roc_curve(class_values,class_labels,0:5);
    [p(ii)] = signrank(class_values(class_labels), class_values(~class_labels)); %#ok
    
    if 0%make figures
        f = figure(...
        'windowstyle', 'normal',...
        'Units', 'pixels',...
        'position', [100 100 500 500],...
        'PaperPositionMode','auto');
        %axes('units', 'pixels', 'position', [0 0 size(mam,2)/2
        %size(mam,1)/2]);
        plot(roc_curve(:,1), roc_curve(:,2), 'LineWidth', 2.0); 
        axis equal; axis([0 1 0 1]);
        set(gca, ...
           'FontName','arial',...
           'FontSize',18,...
           'Xtick', [0 .5 1],...
           'Ytick', [.5 1]);
        text(0.5, 0.4, ['AUC = ', num2str(auc,'%3.2f')], 'FontName','arial', 'FontSize',18);
        roc_name = ['C:\isbe\thesis\figures\11\roc', zerostr(ii,2), '.tif'];
        print('-dtiff', '-noui', '-painters', f, '-r300', roc_name);
        close(f);
    end
    aucs(ii) = auc;
    times(ii) = mean(user_data.ratings(:,3));
    clear user_data;
end
%%
%--------------------------------------------------------------------------
% Compare AUC against experience, job type and time
%
experience = [7 5 2 9 1 4 8 10 6 3]';
rad = [3 5 6 7];
rag = [1 2 4 8 9 10];

display(['Mean AUC = ', num2str(mean(aucs))]);
display(['Mean AUC radiologists = ', num2str(mean(aucs(rad)))]);
display(['Mean AUC = ', num2str(mean(aucs(rag)))]);

[p,h,stats] = ranksum(aucs(rag), aucs(rad));
display(['Radiographers vs radiologists, Mann-whitney U test, p = ', num2str(p)]);
display(stats);
%[p,h,stats] = ranksum(times(rag), times(rad))

[rho,p] = corr(aucs, experience, 'type', 'spearman');
display(['Spearmans correlation between AUC and experienc, rho = ', num2str(rho), ' p = ', num2str(p)]);
[rho,p] = corr(times, experience, 'type', 'spearman');
display(['Spearmans correlation between AUC and time, rho = ', num2str(rho), ' p = ', num2str(p)]);
%[rho,p] = corr(times, aucs, 'type', 'spearman')

%%
%--------------------------------------------------------------------------
% Generate figure of AUC against readers ranked by experience
f = figure(...
        'windowstyle', 'normal',...
        'Units', 'pixels',...
        'position', [100 100 800 500],...
        'PaperPositionMode','auto'); 
bar(experience(rag), aucs(rag), 'b'); hold on; 
bar(experience(rad), aucs(rad), 'r');
set(gca,...
    'Xtick', 1:10, 'XtickLabel', 1:10,...
    'Xlim', [0 11], 'Ylim', [0 1], ...
    'Ygrid', 'on',...
    'FontName','arial',...
    'FontSize',14);
xlabel('Readers ranked by experience', 'FontName','arial', 'FontSize',14);
ylabel('AUC', 'FontName','arial', 'FontSize',14);
title('AUC values for each reader in the observer study', 'FontName','arial', 'FontSize',18);
print('-dtiff', '-noui', '-painters', f, '-r300', 'C:\isbe\thesis\figures\11\experience.tif');
%%
%--------------------------------------------------------------------------
% Analyse user feedback
feedback = zeros(1, 7);
wrong_feedback = zeros(1,7);
f_sum = 0;
f_syn = 0;
f_real = 0;
for ii = 1:10
    display(ii);
    load(['C:\isbe\dev\observer_study\user_data\user', zerostr(ii,2)]);
    f_sum = f_sum + length(user_data.feedback_other);
    for jj = 1:length(user_data.feedback_other)
        if ~isempty(user_data.feedback_other{jj}) && user_data.feedback(jj, 1) <= 30
            display(user_data.feedback_other{jj});
        end
        if user_data.feedback(jj, 1) > 30
            feedback = feedback + user_data.feedback(jj,2:8);
            f_syn = f_syn + 1;
        else
            wrong_feedback = wrong_feedback + user_data.feedback(jj,2:8);
            f_real = f_real + 1;
        end
    end
end
feedback_features = {...
            'No longer think this mass is synthetic',...
            'Mass located incorrectly with respect to breast tissue',...
            'Not enough spicules',...
            'Mass should be more radio-opaque given its size',...
            'Unrealistic texture at mass border',...
            'Mass border too regular',...
            'Mass border too sharply defined'};
        
%%
mass_list = dir('C:\isbe\dev\masses1024x1024\*.mat');
mass_dir = 'C:\isbe\dev\masses1024x1024\';
for ii = 1:30%length(mass_list)
    
    mass = u_load([mass_dir mass_list(ii).name]);
    figure;
    subplot(1,2,1);
    imagesc(mass.mass_ROI); axis image; colormap(gray(256)); colorbar;
    subplot(1,2,2);
    imagesc(mass.background_ROI); axis image; colormap(gray(256)); colorbar;
    
end
%%
syn_10 = u_load('C:\isbe\dev\observer_study\masses\syn_a_mass10.mat');
syn_12 = u_load('C:\isbe\dev\observer_study\masses\syn_a_mass12.mat');
syn_03 = u_load('C:\isbe\dev\observer_study\masses\syn_c_mass03.mat');
syn_06 = u_load('C:\isbe\dev\observer_study\masses\syn_c_mass06.mat');
syn_20 = u_load('C:\isbe\dev\observer_study\masses\syn_c_mass20.mat');
syn_30 = u_load('C:\isbe\dev\observer_study\masses\syn_c_mass30.mat');

write_im_from_colormap(syn_10, 'C:\isbe\thesis\figures\11\syn_10.bmp', gray(256));
write_im_from_colormap(syn_12, 'C:\isbe\thesis\figures\11\syn_12.bmp', gray(256));
write_im_from_colormap(syn_03, 'C:\isbe\thesis\figures\11\syn_03.bmp', gray(256));
write_im_from_colormap(syn_06, 'C:\isbe\thesis\figures\11\syn_06.bmp', gray(256));
write_im_from_colormap(syn_20, 'C:\isbe\thesis\figures\11\syn_20.bmp', gray(256));
write_im_from_colormap(syn_30, 'C:\isbe\thesis\figures\11\syn_30.bmp', gray(256));