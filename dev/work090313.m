% Work Friday 13th March 2009

% I can't believe you're looking again the subtraction method but...

% What if we could do this better by interpolating directly in coarser
% levels of the DT-CWT scaling coefficients of a mass ROI, then reonstructing the estimated background
% using this scaling coefficients and the higher level wavelet coefficients

%%
% First step, look at scaling coefficients in levels 2 to 5 in the mass
% regions

%Load in mass file list
%load C:\isbe\dev\files\u_files.mat

for m = [44 73 78 84 97 99]
    
    mass = u_load(['C:\isbe\dev\masses\', u_files1(m).name]);
    
    [mass_wc mass_sc] = dtwavexfm2(mass.mass_ROI, 6);
    
    figure;
    subplot(2,3,1);
    imagesc(mass.mass_ROI); axis image; colormap(gray(256));
    subplot(2,3,2);
    imagesc(double(mass.mass_ROI) - mass.subtract_ROI); axis image; colormap(gray(256));
    for lev = 2:5
        subplot(2,3,lev+1);
        imagesc(mass_sc{lev}); axis image; colormap(gray(256));
    end
    clear mass;
end

%%
% This new method seems to work quite nicely!!
% More to follow
%load C:\isbe\dev\files\u_files.mat
for m = 141:179
    
    %mass = u_load(['C:\isbe\dev\masses\', u_files1(m).name]);
    mass = u_load(['C:\isbe\dev\annotations\', anno_list(m).name]);
    mass_spline = u_load(['C:\isbe\dev\misc\anno_bg\', anno_list(m).name]);
    [r c] = size(mass.mass_ROI);
    sub_mass = double(mass.mass_ROI) - mass_spline(1:r, 1:c);
    sub_mass(sub_mass < 0) = 0;
%     mass_bw = roipoly(sub_mass, mass.mass_outline(:,1), mass.mass_outline(:,2));
%     for jj = 1:100
%         mass_bw = imdilate(mass_bw, strel('disk', 1));
%     end
%     sub_mass(~mass_bw) = 0;
    
    mass_bg = double(mass.mass_ROI) - sub_mass;
    figure;
    subplot(2,3,1); imagesc(mass.mass_ROI); axis image; colormap(gray(256));
    subplot(2,3,2); imagesc(mass_bg); axis image; colormap(gray(256));
    subplot(2,3,3); imagesc(sub_mass); axis image; colormap(gray(256)); caxis([0 80]);
    hold on;
    plot(mass.mass_outline(:,1), mass.mass_outline(:,2), 'c');
    text(10, 30, num2str(max(sub_mass(:))),'Color',  'r');
    
    
%     mass_spline = u_load(['C:\isbe\dev\misc\anno_bg2\', anno_list(m).name]);
%     [r c] = size(mass.mass_ROI);
%     sub_mass = double(mass.mass_ROI) - mass_spline(1:r, 1:c);
%     mass_bw = roipoly(sub_mass, mass.mass_outline(:,1), mass.mass_outline(:,2));
%     for jj = 1:100
%         mass_bw = imdilate(mass_bw, strel('disk', 1));
%     end
%     sub_mass(~mass_bw) = 0;
    sub_mass = mass.mass_sub_it;
    sub_mass(sub_mass < 0) = 0;
    mass_bg = double(mass.mass_ROI) - double(sub_mass);

    subplot(2,3,4); imagesc(mass.mass_ROI); axis image; colormap(gray(256));
    subplot(2,3,5); imagesc(mass_bg); axis image; colormap(gray(256));
    subplot(2,3,6); imagesc(sub_mass); axis image; colormap(gray(256)); caxis([0 80]);
    hold on;
    plot(mass.mass_outline(:,1), mass.mass_outline(:,2), 'c');
    text(10, 30, num2str(max(sub_mass(:))), 'Color', 'r'); 
    
end
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Let's actually create some models using the new masses
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%First save the new estimated background in the same structure as with the
%original annotated masses
for m = 1:179
    
    %mass = u_load(['C:\isbe\dev\masses\', u_files1(m).name]);
    mass = u_load(['C:\isbe\dev\annotations\', anno_list(m).name]);
    mass_spline = u_load(['C:\isbe\dev\misc\anno_bg\', anno_list(m).name]);
    [r c] = size(mass.mass_ROI);
    sub_mass = double(mass.mass_ROI) - mass_spline(1:r, 1:c);
    sub_mass(sub_mass < 0) = 0;
    mass.mass_sub_it = sub_mass;
    save(['C:\isbe\dev\misc\anno_bg\', anno_list(m).name], 'mass');
end
%%    
% Now turns these into the mass structures used as input to the model:
get_mass_data(anno_list, 'C:\isbe\dev\misc\anno_bg\', 0);

%%
% Now build the model:
load C:\isbe\dev\files\u_files.mat

[mass_model model_id]...
    = generate_mass_AM(u_files1, 'C:\isbe\dev\misc\anno_bg\model_dt_sub_masses',...
    'mass_path', 'C:\isbe\dev\misc\anno_bg\');
%Sort out model weightings

%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


mass_spline7 = subtract_mass_dual_tree('MassList', u_files1(78), 'Level', 6);
mass_spline6 = subtract_mass_dual_tree('MassList', u_files1(78), 'Level', 6, 'n1', 5);
figure; imagesc(double(mass.mass_ROI) - mass_spline7(1:end-1,:)); axis image; colormap(gray(256));
figure; imagesc(mass_spline7); axis image; colormap(gray(256));
mass_spline7o = mass_spline7;
mass_spline7o(mass_spline7o < 0) = 0;
mass_spline6o = mass_spline6;
mass_spline6o(mass_spline6o < 0) = 0;
figure; imagesc(mass_spline6o); axis image; colormap(gray(256));
subtract6 = double(mass.mass_ROI) - mass_spline6(1:end-1,:);
subtract7 = double(mass.mass_ROI) - mass_spline7(1:end-1,:);
subtract7o = subtract7; subtract7o(subtract7o < 0) = 0;
subtract6o = subtract6; subtract7o(subtract6o < 0) = 0;
figure; imagesc(subtract6o); axis image; colormap(gray(256));
subtract7o = subtract7; subtract7o(subtract7o < 0) = 0;
subtract6o = subtract6; subtract6o(subtract6o < 0) = 0;
figure; imagesc(subtract6o); axis image; colormap(gray(256));
figure; imagesc(subtract7o); axis image; colormap(gray(256));
figure; imagesc(double(mass.mass_ROI) - subtract7o(1:end-1,:)); axis image; colormap(gray(256));
figure; imagesc(double(mass.mass_ROI) - subtract7o); axis image; colormap(gray(256));
figure; imagesc(double(mass.mass_ROI) - subtract6o); axis image; colormap(gray(256));
%%