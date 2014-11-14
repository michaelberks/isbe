%Script investigating the use of the STAPLE algorithm to estimate ground
%truth where none is available

%--------------------------------------------------------------------------
%1. First lets recreate the synthetic experiment(s) from the warfield paper

%create 10 BW images for raters each with sensitivity = 0.95 and specificty
% = 0.90
phantom = cell(10,1);
D = zeros(256*256, 10);
figure;
for r = 1:10
    phantom{r} = [rand(256,128) < 0.1 rand(256,128) < 0.95];
    subplot(2,5,r); imagesc(phantom{r}); axis image; colormap gray;
    D(:,r) = phantom{r}(:);
end

%Run STAPLE
[C sens spec it] = staple_consensus(D);
figure; 
subplot(1,2,1); imagesc(reshape(C, 256, 256)); axis image; colormap(gray(256));
subplot(1,2,2); imagesc(reshape(C > .5, 256, 256)); axis image; colormap(gray(256));
display(['Mean sensitivity = ' num2str(mean(sens)) ' +/ ' num2str(std(sens))]);
display(['Mean specificity = ' num2str(mean(spec)) ' +/ ' num2str(std(spec))]);
display(['Convergence reached after ' num2str(it) ' iterations']);
%%
%--------------------------------------------------------------------------
% 2. Now try a similar test but create probabilistic inputs peaked at 0.8
% (for true labels) and 0.2 (for false labels)
phantom = cell(10,1);
D = zeros(256*256, 10);
figure;
for r = 1:10
    phantom{r} = [0.2 + 0.2*randn(256,128) 0.8 + 0.2*randn(256,128)];
    phantom{r}(phantom{r} < 0) = -phantom{r}(phantom{r} < 0);
    phantom{r}(phantom{r} > 1) = 2 - phantom{r}(phantom{r} > 1);
    subplot(2,5,r); imagesc(phantom{r}); axis image; colormap gray;
    D(:,r) = phantom{r}(:);
end

%Run STAPLE
[C sens spec it] = staple_consensus(D);
figure; 
subplot(1,2,1); imagesc(reshape(C, 256, 256)); axis image; colormap(gray(256));
subplot(1,2,2); imagesc(reshape(C > .5, 256, 256)); axis image; colormap(gray(256));
display(['Mean sensitivity = ' num2str(mean(sens)) ' +/ ' num2str(std(sens))]);
display(['Mean specificity = ' num2str(mean(spec)) ' +/ ' num2str(std(spec))]);
display(['Convergence reached after ' num2str(it) ' iterations']);

%Note the much lower sensitivities/specificities of each rater
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%%
%--------------------------------------------------------------------------
% 3. Load in test line images with know ground truth
test_dir = 'M:\chen\data\testimage_contrast0to8_multibars_sin';
prob_dir = 'M:\chen\data\testimage_contrast0to8_multibars_sin\probability_image\';
rater_dirs = dir([prob_dir '18*']);
num_raters = 17;%length(rater_dirs);
D = zeros(512*512, num_raters);

load M:\chen\data\testimage_contrast0to8_multibars_sin\image1.mat
load M:\chen\data\testimage_contrast0to8_multibars_sin\label1.mat
for ii = 1:num_raters
    prob_image = u_load([prob_dir rater_dirs(ii).name '\probability_image001.mat']);
    figure; imagesc(prob_image); axis image; colormap(gray(256));
    D(:,ii) = prob_image(:);
end
%%
[C sens spec it] = staple_consensus(D > .5);
figure; 
subplot(1,2,1); imagesc(reshape(C, 512, 512)); axis image; colormap(gray(256));
subplot(1,2,2); imagesc(reshape(C > .5, 512, 512)); axis image; colormap(gray(256));
display(['Mean sensitivity = ' num2str(mean(sens)) ' +/ ' num2str(std(sens))]);
display(['Mean specificity = ' num2str(mean(spec)) ' +/ ' num2str(std(spec))]);
display(['Convergence reached after ' num2str(it) ' iterations']);
%%
%--------------------------------------------------------------------------
% 4. Now try running the alogithm on some real mammogram inputs
pred{1} = u_load('M:\chen\data\line_detection_mammo\results\predict_mammo_182267.mat');
pred{2} = u_load('M:\chen\data\line_detection_mammo\results\predict_mammo_182269.mat');
pred{3} = u_load('M:\chen\data\line_detection_mammo\results\predict_mammo_182272.mat');
pred{4} = u_load('M:\chen\data\line_detection_mammo\results\predict_mammo_182270.mat');
%
D = cell(5,1);
for ii = 1:5
    figure;
    for jj = 1:4
        subplot(2,2,jj); imagesc(pred{jj}(ii).probability_image); axis image; colormap(gray(256));
        D{ii}(:,jj) = 1-pred{jj}(ii).probability_image(:);
    end
    
    %Run STAPLE
    [C sens spec it] = staple_consensus(D{ii});
    figure; 
    subplot(1,2,1); imagesc(reshape(C, 512, 512)); axis image; colormap(gray(256));
    subplot(1,2,2); imagesc(reshape(C > .5, 512, 512)); axis image; colormap(gray(256));
    display(['Mean sensitivity = ' num2str(mean(sens)) ' +/ ' num2str(std(sens))]);
    display(['Mean specificity = ' num2str(mean(spec)) ' +/ ' num2str(std(spec))]);
    display(['Convergence reached after ' num2str(it) ' iterations']);
end
%%