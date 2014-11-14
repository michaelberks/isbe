%%
%--------------------------------------------------------------------------
%-- Experimental code using snakes to update annotated spicule position
%--------------------------------------------------------------------------
probimage_dir = 'M:\chen\data\predict_masses512x512\';
% probimage_dir =  'M:\chen\data\masses512x512_nms\';
newspicules_dir = 'M:\chen\data\masses512x512_newspicules\';

massimage_dir = 'M:\chen\data\masses512x512\';
spicules_dir = 'M:\chen\data\masses512x512_spicules\';

fg=figure;
frame=1;
Film = struct('cdata', [],'colormap', []);
set(gca,'nextplot','replacechildren');
set(gca,'xlim',[-80 80],'ylim',[-80 80],'NextPlot','replace','Visible','off');
for mass_idx =1:50%179
    %Load in probability images and spicules
    % mass_idx = k; %28
    prob_image = u_load([probimage_dir, 'probability_image', zerostr(mass_idx,3), '.mat']);
    newspicules = u_load([newspicules_dir, 'mass_spicules', zerostr(mass_idx,3), '.mat']);
    spicules = u_load([spicules_dir, 'mass_spicules', zerostr(mass_idx,3), '.mat']);
    mass_image = u_load([massimage_dir, 'mass', zerostr(mass_idx,3), '.mat']);
    feat_img = 1 - prob_image; %Feature image = line probability
    
    clf(fg);
%     figure('windowstyle', 'docked');
    axis off;
    hold on;
    axes('position',[0  0  1/3 1]);
    imagesc(mass_image); axis image; axis off; colormap(gray(256));
    axes('position',[1/3  0  1/3 1]);
    imagesc(mass_image); axis image; axis off; colormap(gray(256));
    axes('position',[2/3  0  1/3 1]);
    imagesc(feat_img); axis image; axis off; colormap(gray(256));
    truesize;
    pause(0.5);
    Film(frame) = getframe(gcf);
    frame=frame+1;
    
    clf(fg);
% figure('windowstyle', 'docked');
    axis off;
    hold on;
    axes('position',[0  0  1/3 1]);
    imagesc(mass_image); axis image; axis off; colormap(gray(256));
    hold on;
    if ~isempty(spicules)
        for ii = 1:length(spicules)
            
            %Select i-th spicule
            spicule = spicules{ii};
            
            plot(spicule(:,1), spicule(:,2),  'g:', 'LineWidth',1);
        end
    end
    
    
    axes('position',[1/3  0  1/3 1]);
    imagesc(mass_image); axis image; axis off; colormap(gray(256));
    hold on;
    if ~isempty(newspicules)
        for ii = 1:length(newspicules)
            
            %Select i-th spicule
            spicule = newspicules{ii};
            
            plot(spicule(:,1), spicule(:,2),  'r:', 'LineWidth',1);
        end
    end
    axes('position',[2/3  0  1/3 1]);
    imagesc(feat_img); axis image; axis off; colormap(gray(256));
    hold on;
    if ~isempty(newspicules)
        for ii = 1:length(newspicules)
            
            %plot original spicule
            spicule = spicules{ii};
            
            plot(spicule(:,1), spicule(:,2),  'g:', 'LineWidth',1);
            
            %Select i-th spicule
            spicule = newspicules{ii};
            
            plot(spicule(:,1), spicule(:,2),  'r:', 'LineWidth',1);
        end
    end
    %             truesize
    Film(frame) = getframe(gcf);
    pause(0.5);
    frame=frame+1;
    
    disp(mass_idx);
    clear prob_image newspicules spicules mass_image;
end

%% create a vedio of the results
% create_avi(Film, [newspicules_dir, 'spicule_snakes1_50.avi'], 5);

