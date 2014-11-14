%%
%--------------------------------------------------------------------------
%-- Experimental code using snakes to update annotated spicule position
%--------------------------------------------------------------------------

n = 102;
bin_edges = linspace(0,1,n);

%Load in probability images and spicules and generate sample of data
line_strength_probs = []; %We're going to grow this in a loop, hopefully the inefficiency won't hurt too much...
spicule_sample_pts = cell(179, 1);
spic_pts_per_image = zeros(179,1);
mass_pdf = 0;
num_spics = 0; 
for mass_idx = 1:179

    prob_image = 1-u_load(['M:\chen\data\predict_masses512x512\probability_image', zerostr(mass_idx,3), '.mat']);
    spicules = u_load(['M:\chen\data\masses512x512_newspicules\mass_spicules', zerostr(mass_idx,3), '.mat']);

    [row col] = size(prob_image); %should always be 512x512
    
    mass_pdf = mass_pdf + histc(prob_image(:), bin_edges);
    num_spics = num_spics + length(spicules);
    
    %figure; imagesc(prob_image); axis image; colormap(gray(256)); hold on;
    for ii = 1:length(spicules)
        
        %get i-th spicule
        spic = spicules{ii};
        
        %throwaway any points outside of image
        discard_pts = any([spic < 1, spic(:,1) > col, spic(:,2) > col],2);
        spic(discard_pts,:) = [];
        
        %plot(spic(:,1), spic(:,2), 'g');
        
        spic_idx = sub2ind([row col], spic(:,2), spic(:,1));
        line_strength_probs = [line_strength_probs; prob_image(spic_idx(:))]; %#ok
        spicule_sample_pts{mass_idx} = [spicule_sample_pts{mass_idx}; spic];
    end
    spic_pts_per_image(mass_idx) = size(spicule_sample_pts{mass_idx}, 1);
end
%
spicule_pdf = histc(line_strength_probs, bin_edges);
figure; bar(bin_edges, spicule_pdf); 
title('Distribution of line probabilities in mass spicules');
figure; bar(bin_edges, mass_pdf);
title('Distribution of line probabilities in all mass regions');

spic_pts_cumsum = cumsum(spic_pts_per_image);
%
%--------------------------------------------------------------------------
% Generate sample of points to take from normal data
%%
%First loop through the normal data and get a histogram of line probs
norm_pdf = 0;
for norm_idx = 1:89

    prob_image = 1-u_load(['M:\chen\data\predict_normal_512x512\probability_image', zerostr(norm_idx,3), '.mat']);
    norm_pdf = norm_pdf + histc(prob_image(:), bin_edges);
end
figure; bar(bin_edges, norm_pdf);
title('Distribution of line probabilities in all normal regions');
%%
%So now we know, for the j-th n bins, there are in total norm_pdf(j) pts
% of which we want to select spicule_pdf(j).
bin_idx = cell(n,1);
bin_totals = zeros(n,1);
for jj = 1:n
    %get random permutation of norm_pdf(j)
    rp = randperm(norm_pdf(jj));
    
    %select the first spicule_pdf(jj) and sort
    bin_idx{jj} = sort(rp(1:spicule_pdf(jj)));
    bin_totals(jj) = length(bin_idx{jj});
end
display(sum(bin_totals));
%%
%Now go through each normal image and for each bin, work out which points
%fall in the bin and then whether we have chosen to select them
norm_line_probs = [];
normal_sample_pts = cell(89, 1);

norm_pts_per_image = zeros(89, 1);
bin_totals2 = zeros(n, 1);

for norm_idx = 1:89

    prob_image = 1-u_load(['M:\chen\data\predict_normal_512x512\probability_image', zerostr(norm_idx,3), '.mat']);
    
    %for each bin
    for jj = 1:n
        
        if jj < n
            %find pts in this image that fall in bin
            possible_idx = find((prob_image >= bin_edges(jj)) & (prob_image < bin_edges(jj+1)));
        else
            possible_idx = find(bin_edges(jj) == prob_image);
        end

        %work out which pre-selected samples apply to this image
        no_pts_in_image = sum(bin_idx{jj} <= length(possible_idx));
            
        if no_pts_in_image
            
            bin_image_idx = possible_idx(bin_idx{jj}(1:no_pts_in_image));

            norm_pts_per_image(norm_idx) = norm_pts_per_image(norm_idx) + no_pts_in_image;
            bin_totals2(jj) = bin_totals2(jj) + no_pts_in_image;
            
            %remove these points from bin_idx
            bin_idx{jj}(1:no_pts_in_image) = [];
            

            %Get sample so we can check at the end we've done this right
            norm_line_probs = [norm_line_probs; prob_image(bin_image_idx(:))]; %#ok

            %Save the sample points
            [line_pts_y line_pts_x] = ind2sub(size(prob_image), bin_image_idx(:));
            normal_sample_pts{norm_idx} = [normal_sample_pts{norm_idx}; line_pts_x line_pts_y];
        end
        
        %subtract so counter ready for next image
        bin_idx{jj} = bin_idx{jj} - length(possible_idx);
    end

end

norm_pts_cumsum = cumsum(norm_pts_per_image);
display(sum(norm_pts_per_image));
display(sum(bin_totals2));

%Check the sample of line probs looks like we thought it should
norm_spicule_pdf = histc(norm_line_probs, bin_edges);
figure; bar(bin_edges, norm_spicule_pdf);
title('Distribution of line probabilities in selected normal points');
hold on;
plot(bin_edges, bin_totals, 'g');
plot(bin_edges, bin_totals2, 'r:');

%%
mkdir C:\isbe\dev\classification\data\spicules
save C:\isbe\dev\classification\data\spicules\experiment_data.mat normal_sample_pts spicule_sample_pts bin_idx bin_edges *pdf *pts_cumsum
%%
load C:\isbe\dev\classification\data\spicules\experiment_data.mat
%%
%Check the spicule points lie on spicules
for ii = 151:179
    mass = u_load(['C:\isbe\dev\image_data\masses512x512\mass', zerostr(ii,3), '.mat']);
    spicules = u_load(['C:\isbe\dev\image_data\masses512x512_spicules\mass_spicules', zerostr(ii,3), '.mat']);
    
    if ~isempty(spicule_sample_pts{ii})
        figure; imagesc(mass); axis image; colormap(gray(256)); hold on;
        for jj = 1:length(spicules)
            plot(spicules{jj}(:,1), spicules{jj}(:,2), 'g');
        end
        plot(spicule_sample_pts{ii}(:,1), spicule_sample_pts{ii}(:,2), 'r.');
    end
end
%%
%Masses are ok, now have  a look where the points in normal images lie        
for ii = 61:89
    norm = u_load(['C:\isbe\dev\image_data\normal_512\bg', zerostr(ii,3), '.mat']);
    
    if ~isempty(normal_sample_pts{ii})
        figure; imagesc(norm); axis image; colormap(gray(256)); hold on;
        plot(normal_sample_pts{ii}(:,1), normal_sample_pts{ii}(:,2), 'r.');
    end
end
%%
%--------------------------------------------------------------------------
% Sample different set of data for both the spicule and normal sample
% points
load C:\isbe\dev\classification\data\spicules\experiment_data.mat

mass_list = dir('C:\isbe\dev\image_data\masses512x512\*.mat');
norm_list = dir('C:\isbe\dev\image_data\normal_512\*.mat');

mass_data = sample_spicule_training_data(...
    'image_dir', 'C:\isbe\dev\image_data\masses512x512\',...
    'image_sample_pts', spicule_sample_pts,...
    'image_list', mass_list,...
    'plot', 0);

norm_data = sample_spicule_training_data(...
    'image_dir', 'C:\isbe\dev\image_data\normal_512\',...
    'image_sample_pts', normal_sample_pts,...
    'image_list', norm_list,...
    'plot', 0);

save C:\isbe\dev\classification\data\spicules\experiment_samples.mat mass_data norm_data
clear; pack;
%%
load C:\isbe\dev\classification\data\spicules\experiment_samples.mat
random_forest = mb_random_forest_class_train_boot(...
    'X', [mass_data; norm_data],...
    'y', [true(size(mass_data,1),1); true(size(norm_data,1),1)],...
    'n_trees', 2,...
    'do_oob', 0);    

%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Cross fold forests have been run on the cluster, load in the results
figure; hold all;
%colors = 'xxrgbk';
leg_text = cell(1,1);
for level = 3
    pooled_test_labels = [];
    pooled_test_scores = [];
    fold_aucs = zeros(10,1);
%     figure;
    for ii = [1 4:10]
        load(['C:\isbe\dev\classification\data\spicules\dt\rf_spic_dt_3_', zerostr(level,1), '_all_', zerostr(ii,2), '_results.mat']);
        fold_accuracy = mean(test_labels == str2double(test_predictions));
        display(['accuracy of fold ', num2str(ii), ' = ', num2str(fold_accuracy)]);

        test_scores = test_votes(:,2) / sum(test_votes(1,:));

        pooled_test_labels = [pooled_test_labels; test_labels]; %#ok
        pooled_test_scores = [pooled_test_scores; test_scores]; %#ok

%         [roc_pts fold_aucs(ii)] = calculate_roc_curve(test_scores,test_labels,linspace(-0.0001,1.0001,100));
%         subplot(2,5,ii); plot(roc_pts(:,1), roc_pts(:,2));
%         xlabel('False positive rate');
%         ylabel('True positive rate');
%         title(['ROC curve: ', num2str(ii)]);
%         text(0.5, 0.4, ['AUC = ', num2str(fold_aucs(ii))]);
    end 
    %
    [pooled_roc_pts pooled_auc] = calculate_roc_curve(pooled_test_scores,pooled_test_labels,linspace(-0.0001,1.0001,100));


%     figure;
    %plot(pooled_roc_pts(:,1), pooled_roc_pts(:,2), colors(level));
    plot(pooled_roc_pts(:,1), pooled_roc_pts(:,2));
    xlabel('False positive rate');
    ylabel('True positive rate');
    title('ROC curve for pooled cross-fold classification');
    text(0.5, 0.1 + level/10, ['AUC = ', num2str(pooled_auc)]);
    
    leg_text{level-2} = [num2str(level) ' levels, A_{z}=' num2str(pooled_auc,3)];
end
legend(leg_text, 'location', 'southeast');
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%%
save C:\isbe\dev\classification\data\spicules\experiment_results.mat fold_aucs pooled_*
%%
%--------------------------------------------------------------------------
%Now look at how things are classified on real images
mass_idx = 28;
norm_idx = 6;
mass = u_load(['C:\isbe\dev\image_data\masses512x512\mass', zerostr(mass_idx,3), '.mat']);
norm = u_load(['C:\isbe\dev\image_data\normal_512\norm', zerostr(norm_idx,3), '.mat']);
load C:\isbe\dev\classification\data\spicules\rf_all_data_3_6_all_10.mat
[mass_prob_spic] = classify_image(...
    'image_in', mass, ...
     'forest', random_forest,...
     'forest_type', 'isbe_boot');
 
[norm_prob_spic] = classify_image(...
    'image_in', norm, ...
     'forest', random_forest,...
     'forest_type', 'isbe_boot');
 
mass_prob_line = 1-u_load(['M:\chen\data\predict_masses512x512\probability_image', zerostr(mass_idx,3), '.mat']);
norm_prob_line = 1-u_load(['M:\chen\data\predict_normal_512x512\probability_image', zerostr(norm_idx,3), '.mat']);

mass_combined = mass_prob_line .* exp(i*mass_prob_spic*pi);
norm_combined = norm_prob_line .* exp(i*norm_prob_spic*pi);
figure;
subplot(1,2,1); image(complex2rgb(mass_combined)); axis image;
subplot(1,2,2); image(complex2rgb(norm_combined)); axis image;
%%
% Run the rest of the classifications on the cluster, now we can look at
% any by...
for mass_idx = 1:179
    try
        mass = u_load(['C:\isbe\dev\image_data\masses512x512\mass', zerostr(mass_idx,3), '.mat']);
        mass_prob_line = 1-u_load(['M:\chen\data\predict_masses512x512\probability_image', zerostr(mass_idx,3), '.mat']);
        mass_prob_spic = 1-u_load(['C:\isbe\dev\classification/data/masses_512_spic_probs_6/prob_spic', zerostr(mass_idx,3), '.mat']);
        mass_combined = complex2rgb(mass_prob_line .* exp(i*mass_prob_spic*pi));
        mass_prob = complex2rgb(exp(i*mass_prob_spic*pi));
        mass_thresh = complex2rgb(exp(i*(mass_prob_spic > .5)*pi));
        
        orig_rgb = cat(3, mass, mass, mass) / 255;
 
        imwrite(cat(2, orig_rgb, mass_combined), ['M:\asymmetry_project\experiments\spicule_classification\figures\masses\spic_prob', zerostr(mass_idx,3), '_a.bmp']);
        imwrite(cat(2, orig_rgb, mass_prob), ['M:\asymmetry_project\experiments\spicule_classification\figures\masses\spic_prob', zerostr(mass_idx,3), '_b.bmp']);
        imwrite(cat(2, orig_rgb, mass_thresh), ['M:\asymmetry_project\experiments\spicule_classification\figures\masses\spic_prob', zerostr(mass_idx,3), '_c.bmp']);
%         figure; 
%         subplot(1,2,1); image(complex2rgb(mass_combined)); axis image;
%         subplot(1,2,2); image(complex2rgb(mass_thresh)); axis image;

        %imwrite(complex2rgb([mass_combined mass_thresh]), ['C:\isbe\dev\classification\data\spicules\figures\masses\spic_prob', zerostr(mass_idx,3), '.bmp']);
        %imwrite(complex2rgb([exp(i*mass_prob_spic*pi) exp(i*(mass_prob_spic > .5)*pi)]), ['C:\isbe\dev\classification\data\spicules\figures\masses\spic_prob_unmod', zerostr(mass_idx,3), '.bmp']);
        %imwrite(complex2rgb(mass_thresh), ['C:\isbe\dev\classification\data\spicules\figures\masses\spic_class', zerostr(mass_idx,3), '.bmp']);
    catch
        display(['image ', num2str(mass_idx), ' not found']);
    end
end
%%
for norm_idx = 1:89
    try
        norm = u_load(['C:\isbe\dev\image_data\normal_512\norm', zerostr(norm_idx,3), '.mat']);
        norm_prob_line = 1-u_load(['M:\chen\data\predict_normal_512x512\probability_image', zerostr(norm_idx,3), '.mat']);
        norm_prob_spic = 1-u_load(['C:\isbe\dev\classification/data/normal_512_spic_probs_6/prob_spic', zerostr(norm_idx,3), '.mat']);
        norm_combined = complex2rgb(norm_prob_line .* exp(i*norm_prob_spic*pi));
        norm_prob = complex2rgb(exp(i*norm_prob_spic*pi));
        norm_thresh = complex2rgb(exp(i*(norm_prob_spic > .5)*pi));
        
        orig_rgb = cat(3, norm, norm, norm) / 255;
 
        imwrite(cat(2, orig_rgb, norm_combined), ['M:\asymmetry_project\experiments\spicule_classification\figures\normals\spic_prob', zerostr(norm_idx,3), '_a.bmp']);
        imwrite(cat(2, orig_rgb, norm_prob), ['M:\asymmetry_project\experiments\spicule_classification\figures\normals\spic_prob', zerostr(norm_idx,3), '_b.bmp']);
        imwrite(cat(2, orig_rgb, norm_thresh), ['M:\asymmetry_project\experiments\spicule_classification\figures\normals\spic_prob', zerostr(norm_idx,3), '_c.bmp']);
%         figure; 
%         subplot(1,2,1); image(complex2rgb(norm_combined)); axis image;
%         subplot(1,2,2); image(complex2rgb(norm_thresh)); axis image;
        %imwrite(complex2rgb([norm_combined norm_thresh]), ['C:\isbe\dev\classification\data\spicules\figures\normals\spic_prob', zerostr(norm_idx,3), '.bmp']);
        %imwrite(complex2rgb(norm_thresh), ['C:\isbe\dev\classification\data\spicules\figures\normals\spic_class', zerostr(norm_idx,3), '.bmp']);
        
    catch
        display(['image ', num2str(norm_idx), ' not found']);
    end
end
%%

for fold = 1:10
    if fold == 1
        start_m = 1;
        start_m_i = 1;
    else
        start_m = find((spic_pts_cumsum / spic_pts_cumsum(end))>= (fold-1)/n_fold, 1) + 1; %#ok
        start_m_i = spic_pts_cumsum(start_m-1) + 1;
    end
    end_m = find((spic_pts_cumsum / spic_pts_cumsum(end)) >= fold/n_fold, 1); %#ok
    end_m_i = spic_pts_cumsum(end_m);
    display([start_m end_m start_m_i end_m_i]);
    
%     if fold == 1
%         start_n = 1;
%         start_n_i = 1;
%     else
%         start_n = find((norm_pts_cumsum / norm_pts_cumsum(end))>= (fold-1)/n_fold, 1) + 1; %#ok
%         start_n_i = norm_pts_cumsum(start_n-1) + 1;
%     end
%     end_n = find((norm_pts_cumsum / norm_pts_cumsum(end)) >= fold/n_fold, 1); %#ok
%     end_n_i = norm_pts_cumsum(end_n); %#ok
%     display([start_n end_n start_n_i end_n_i]);
end
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Figure for IWDM paper
norm_idx = 68;

norm = u_load(['C:\isbe\dev\image_data\normal_512\norm', zerostr(norm_idx,3), '.mat']);
norm_prob_line = 1-u_load(['M:\chen\data\predict_normal_512x512\probability_image', zerostr(norm_idx,3), '.mat']);
norm_prob_spic = 1-u_load(['C:\isbe\dev\classification/data/normal_512_spic_probs_6/prob_spic', zerostr(norm_idx,3), '.mat']);

norm_combined = complex2rgb(norm_prob_line .* exp(i*norm_prob_spic*pi));
norm_prob = complex2rgb(exp(i*norm_prob_spic*pi));
norm_thresh = complex2rgb(exp(i*(norm_prob_spic > .5)*pi));

write_im_from_colormap(norm,...
    ['K:\isbe\conferences_and_symposia\iwdm2010\bars\figures\norm', zerostr(norm_idx,3), '.bmp'],...
    gray(256));

write_im_from_colormap(norm_prob_line,...
    ['K:\isbe\conferences_and_symposia\iwdm2010\bars\figures\line_prob_norm', zerostr(norm_idx,3), '.bmp'],...
    gray(256), [0 1]);

imwrite(norm_combined,...
    ['K:\isbe\conferences_and_symposia\iwdm2010\bars\figures\spic_prob_norm', zerostr(norm_idx,3), '_a.bmp']);
imwrite(norm_prob,...
    ['K:\isbe\conferences_and_symposia\iwdm2010\bars\figures\spic_prob_norm', zerostr(norm_idx,3), '_b.bmp']);
imwrite(norm_thresh,...
    ['K:\isbe\conferences_and_symposia\iwdm2010\bars\figures\spic_prob_norm', zerostr(norm_idx,3), '_c.bmp']);
%%
mass_idx = 46;

mass = u_load(['C:\isbe\dev\image_data\masses512x512\mass', zerostr(mass_idx,3), '.mat']);
mass_prob_line = 1-u_load(['M:\chen\data\predict_masses512x512\probability_image', zerostr(mass_idx,3), '.mat']);
mass_prob_spic = 1-u_load(['C:\isbe\dev\classification/data/masses_512_spic_probs_6/prob_spic', zerostr(mass_idx,3), '.mat']);

mass_combined = complex2rgb(mass_prob_line .* exp(i*mass_prob_spic*pi));
mass_prob = complex2rgb(exp(i*mass_prob_spic*pi));
mass_thresh = complex2rgb(exp(i*(mass_prob_spic > .5)*pi));

write_im_from_colormap(mass,...
    ['K:\isbe\conferences_and_symposia\iwdm2010\bars\figures\mass', zerostr(mass_idx,3), '.bmp'],...
    gray(256));

write_im_from_colormap(mass_prob_line,...
    ['K:\isbe\conferences_and_symposia\iwdm2010\bars\figures\line_prob_mass', zerostr(mass_idx,3), '.bmp'],...
    gray(256), [0 1]);

imwrite(mass_combined,...
    ['K:\isbe\conferences_and_symposia\iwdm2010\bars\figures\spic_prob_mass', zerostr(mass_idx,3), '_a.bmp']);
imwrite(mass_prob,...
    ['K:\isbe\conferences_and_symposia\iwdm2010\bars\figures\spic_prob_mass', zerostr(mass_idx,3), '_b.bmp']);
imwrite(mass_thresh,...
    ['K:\isbe\conferences_and_symposia\iwdm2010\bars\figures\spic_prob_mass', zerostr(mass_idx,3), '_c.bmp']);