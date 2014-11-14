anno_list = dir('C:\isbe\dev\annotations\*.mat');
mass_list_1024 = dir('C:\isbe\dev\masses1024x1024\*.mat');
mass_list_512 = dir('C:\isbe\dev\image_data\masses512x512\*.mat');
mass_list = dir('C:\isbe\dev\masses\*.mat');
%
n_o_s = 0;
for ii = [28 47]%length(anno_list)
    anno = u_load(['C:\isbe\dev\annotations\', anno_list(ii).name]);
    mass1024 = u_load(['C:\isbe\dev\masses1024x1024\', mass_list_1024(ii).name]);
    mass512 = u_load(['C:\isbe\dev\image_data\masses512x512\', mass_list_512(ii).name]);
    mass = u_load(['C:\isbe\dev\masses\', mass_list(ii).name]);
    
    if ~isempty(anno.mass_spicules)
        %n_o_s = n_o_s + length(mass.mass_spicules);
        figure('name', mass_list(ii).name); 
        subplot(1,2,1); imagesc(anno.mass_ROI); axis image; colormap(gray(256)); hold on;
        plot(anno.mass_outline(:,1), anno.mass_outline(:,2));
        for jj = 1:length(anno.mass_spicules)
            spic = anno.mass_spicules(jj).outline;
            plot(spic(:,1),spic(:,2), 'r--');
        end
        
        subplot(1,2,2); imagesc(mass512); axis image; colormap(gray(256)); hold on;
        for jj = 1:length(mass1024.mass_spicules)
            spic = mass1024.mass_spicules(jj).outline;
            spic = ((spic - 1)/2)+1;
            plot(spic(:,1),spic(:,2), 'r--');
        end
        
    end
    clear mass;
end
%%
mkdir C:\isbe\dev\image_data\masses512x512_spicules\;
mass_list_1024 = dir('C:\isbe\dev\masses1024x1024\*.mat');
for ii = 1:length(mass_list_1024);
    mass1024 = u_load(['C:\isbe\dev\masses1024x1024\', mass_list_1024(ii).name]);
    spicules = cell(length(mass1024.mass_spicules),1);
    for jj = 1:length(mass1024.mass_spicules)
        spic = mass1024.mass_spicules(jj).outline;
        spic = ((spic - 1)/2)+1;
        dists = cumsum([0; sqrt(sum(diff(spic).^2, 2))]);
        spic = interp1(dists, spic, 0:10:dists(end), 'linear');
        
        %Throw away points off the edge of the image
        discard_idx = any(spic < 1, 2) | any(spic > 512, 2);
        spic(discard_idx,:) = [];
        spicules{jj} = spic;
    end
    save(['C:\isbe\dev\image_data\masses512x512_spicules\mass_spicules', zerostr(ii,3), '.mat'], 'spicules');
    clear mass1024 spicules
end
%%
%%
mass_roi1 = u_load('C:\isbe\dev\image_data\masses512x512\mass028.mat');
mass_roi2 = u_load('C:\isbe\dev\image_data\masses512x512\mass047.mat');
mass_roi1_prob = u_load('C:\isbe\dev\image_data\line_detection_mammo\mass_roi1_prob.mat');
mass_roi2_prob = u_load('C:\isbe\dev\image_data\line_detection_mammo\mass_roi2_prob.mat');

mass_roi1_spicules = u_load('C:\isbe\dev\image_data\masses512x512_spicules\mass_spicules028.mat');
mass_roi2_spicules = u_load('C:\isbe\dev\image_data\masses512x512_spicules\mass_spicules047.mat');
%
figure;
subplot(1,2,1); imagesc(mass_roi1); axis image; colormap(gray(256)); hold on;
for jj = 1:length(mass_roi1_spicules)
    spic = mass_roi1_spicules{jj};
    plot(spic(:,1),spic(:,2), 'r--');
    plot(spic(:,1),spic(:,2), 'bx');
end
subplot(1,2,2); imagesc(1-mass_roi1_prob); axis image; colormap(gray(256)); hold on;
for jj = 1:length(mass_roi1_spicules)
    spic = mass_roi1_spicules{jj};
    plot(spic(:,1),spic(:,2), 'r--');
    plot(spic(:,1),spic(:,2), 'bx');
end
%
figure;
subplot(1,2,1); imagesc(mass_roi2); axis image; colormap(gray(256)); hold on;
for jj = 1:length(mass_roi1_spicules)
    spic = mass_roi2_spicules{jj};
    plot(spic(:,1),spic(:,2), 'r--');
    plot(spic(:,1),spic(:,2), 'bx');
end
subplot(1,2,2); imagesc(1-mass_roi2_prob); axis image; colormap(gray(256)); hold on;
for jj = 1:length(mass_roi2_spicules)
    spic = mass_roi2_spicules{jj};
    plot(spic(:,1),spic(:,2), 'r--');
    plot(spic(:,1),spic(:,2), 'bx');
end
%%
%--------------------------------------------------------------------------
%-- Experimental code using snakes to update annotated spicule position
%--------------------------------------------------------------------------

%Load in probability images and spicules
ii = 28;
prob_image = u_load(['M:\chen\data\predict_masses512x512\probability_image', zerostr(ii,3), '.mat']);
spicules = u_load(['M:\chen\data\masses512x512_spicules\mass_spicules', zerostr(ii,3), '.mat']);

%Set snake parameters
alpha = 0; %Weight to minimise length (set to zero as don't want to shrink spicule)
beta = 0.1; %Weight to minimise curvature
max_delta_x = 3; %Num pixels range horizontally
resol_x = 1; 
max_delta_y = 3; %Num pixels range vertically
resol_y = 1; 
feat_img = 1 - prob_image; %Feature image = line probability

%
for ii = 1:length(spicules)
    spicule = spicules{ii};
    e_old = Inf;
    e_new = 1e6;
    thresh = 1e-6;
    iter = 1;
    while e_new+thresh < e_old;
        figure; imagesc(feat_img); axis image; colormap(gray(256)); hold on;
        plot(spicule(:,1), spicule(:,2), 'g');

        e_old = e_new;
        [spicule, e_new] = mb_snake(spicule, alpha, beta, max_delta_x, resol_x, max_delta_y, resol_y, feat_img);

        plot(spicule(:,1), spicule(:,2), 'bx');
        plot(spicule(:,1), spicule(:,2), 'r');

        display(['Iteration ', num2str(iter), ': energy = ', num2str(e_new)]);
        iter = iter + 1;
    end
end