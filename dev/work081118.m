% Let's test to see whether we can simply throw away all the information
% from non-maximal sub-bands

% Let's look at Lenna

lenna = u_load('C:\isbe\matlab_code\trunk\wavelet_dtcwt_toolbox4_1\lenna.mat');

dt_lenna = dtwavexfm2(lenna, 4, 'near_sym_b','qshift_b');

%%
% Try throwing away all data from non-maximal subbands
dt_lenna_max = dt_lenna;
for lev = 1:4
    
    [dummy, max_subband_idx] = max(dt_lenna_max{lev}, [], 3);
    
    for subband = 1:6
        %for each subband, set to zero all coeffcients that aren't maximal
        %in this subband
        
        temp = dt_lenna_max{lev}(:,:,subband);
        temp(max_subband_idx ~= subband) = 0;
        
        dt_lenna_max{lev}(:,:,subband) = temp;
        
    end
end

% Now reconstruct the image
lenna_max = dtwaveifm2(dt_lenna_max, 'near_sym_b','qshift_b');

%%
% Try keeping data from maximal and adjacent sub-bands
dt_lenna_max2 = dt_lenna;
for lev = 1:3
    
    [dummy, max_subband_idx] = max(dt_lenna_max{lev}, [], 3);
    
    for subband = 1:6
        %for each subband, set to zero all coeffcients that aren't maximal
        %in this subband
        if subband == 6
            sub_up = 1;
        else
            sub_up = subband+1;
        end
        if subband == 1
            sub_down = 6;
        else
            sub_down = subband-1;
        end
        temp = dt_lenna_max2{lev}(:,:,subband);
        temp(max_subband_idx ~= subband & max_subband_idx ~= sub_up & max_subband_idx ~= sub_down) = 0;
        
        dt_lenna_max2{lev}(:,:,subband) = temp;
        
    end
end

% Now reconstruct the image
lenna_max2 = dtwaveifm2(dt_lenna_max2, 'near_sym_b','qshift_b');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lets look at how the magnitudes of points vary across sub-bands
% Do level 4 first as not much data

for lev = 1:4
    figure;
    hold on;
    title(['Line plot of magnitude across sub-bands - Level ', num2str(lev)]);

    [dummy, max_subband_idx] = max(dt_lenna{lev}, [], 3);
    dt_lev = reshape(dt_lenna{lev}, [], 6);
    for subband = 1:6
        temp_idx = reshape(max_subband_idx == subband, [], 1);
        temp = circshift(dt_lev(temp_idx, :), [0 3 - subband]);
        %subplot(2,3,subband); 
        plot(abs(temp)');
    end
end
%%
dt_lenna2 = cell(4,1);
for lev = 1:4

    [dummy, max_subband_idx] = max(dt_lenna{lev}, [], 3);
    dt_lenna2{lev} = reshape(dt_lenna{lev}, [], 6);
    for subband = 1:6
        temp_idx = reshape(max_subband_idx == subband, [], 1);
        dt_lenna2{lev}(temp_idx, :) = circshift(dt_lenna2{lev}(temp_idx, :), [0 1 - subband]);
        %subplot(2,3,subband); 
        
    end
    
end
%%
for lev = 1:4
    figure;
    for col = 1:6
        dt_lev_abs = abs(dt_lenna2{lev}(:,col));
        subplot(2,3,col); hist(dt_lev_abs);
        title(['DT magnitudes - Level ', num2str(lev), ' Sub-band ', num2str(col)]);
    end
    figure;
    for col = 1:6
        dt_lev_abs_percent = abs(dt_lenna2{lev}(:,col)) ./ abs(dt_lenna2{lev}(:,1));
        subplot(2,3,col); hist(dt_lev_abs_percent);
        title(['DT percentage of maximum magnitudes - Level ', num2str(lev), ' Sub-band ', num2str(col)]);
    end
end

for lev = 1:4
    figure;
    for col = 1:6
        dt_lev_angle = angle(dt_lenna2{lev}(:,col));
        subplot(2,3,col); hist(dt_lev_angle);
        title(['DT angles - Level ', num2str(lev), ' Sub-band ', num2str(col)]);
    end
    figure;
    for col = 1:6
        dt_lev_angle_diff = angle(dt_lenna2{lev}(:,col) .* conj(dt_lenna2{lev}(:,1)));
        subplot(2,3,col); hist(dt_lev_angle_diff);
        title(['DT difference angles from maximum magnitudes - Level ', num2str(lev), ' Sub-band ', num2str(col)]);
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Try reconstructing lenna with original magnitudes but random pahse
%%
% Try throwing away all data from non-maximal subbands
dt_lenna_rp = dt_lenna;
for lev = 1:4
    
    % 1) check we're assigning correctly below
    %dt_lenna_rp{lev} = abs(dt_lenna{lev}) .* exp(i*angle(dt_lenna{lev}));
    
    % 2) try zero phase
    %dt_lenna_rp{lev} = abs(dt_lenna{lev});
    
    % 3) try random phase
    dt_lenna_rp{lev} = abs(dt_lenna{lev}) .* exp(i*2*pi*rand(size(dt_lenna{lev})));
 
    [dummy, max_subband_idx] = max(dt_lenna{lev}, [], 3);
    
    for subband = 1:6
        %for each subband, set to zero all coeffcients that aren't maximal
        %in this subband
        
        temp = dt_lenna{lev}(:,:,subband);
        temp_rp = dt_lenna_rp{lev}(:,:,subband);
        temp(max_subband_idx ~= subband) = temp_rp(max_subband_idx ~= subband);
        
        dt_lenna_rp{lev}(:,:,subband) = temp;
        
    end
end

% Now reconstruct the image
lenna_rp = dtwaveifm2(dt_lenna_rp, 'near_sym_b','qshift_b');
figure; imagesc(lenna_rp); axis image; colormap(gray(256));

%%
% Try keeping data from maximal and adjacent sub-bands
dt_lenna_rp2 = dt_lenna;
for lev = 1:4
    
    % 1) check we're assigning correctly below
    %dt_lenna_rp{lev} = abs(dt_lenna{lev}) .* exp(i*angle(dt_lenna{lev}));
    
    % 2) try zero phase
    dt_lenna_rp2{lev} = abs(dt_lenna{lev});
    
    % 3) try random phase
    %dt_lenna_rp2{lev} = abs(dt_lenna{lev}) .* exp(i*2*pi*rand(size(dt_lenna{lev})));
 
    [dummy, max_subband_idx] = max(dt_lenna{lev}, [], 3);
    
    for subband = 1:6
        %for each subband, set to zero all coeffcients that aren't maximal
        %in this subband
        
        %for each subband, set to zero all coeffcients that aren't maximal
        %in this subband
        if subband == 6
            sub_up = 1;
        else
            sub_up = subband+1;
        end
        if subband == 1
            sub_down = 6;
        else
            sub_down = subband-1;
        end
        
        temp = dt_lenna{lev}(:,:,subband);
        temp_rp = dt_lenna_rp2{lev}(:,:,subband);
        
        temp_idx = max_subband_idx ~= subband & max_subband_idx ~= sub_up & max_subband_idx ~= sub_down;
        temp(temp_idx) = temp_rp(temp_idx);
        
        dt_lenna_rp2{lev}(:,:,subband) = temp;
        
    end
end

% Now reconstruct the image
lenna_rp2 = dtwaveifm2(dt_lenna_rp2, 'near_sym_b','qshift_b');
figure; imagesc(lenna_rp2); axis image; colormap(gray(256));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Try reconstructing lenna with original phases but random magnitudes
%%
% Try throwing away all data from non-maximal subbands
dt_lenna_rm = dt_lenna;
for lev = 1:4
    
    [max_lenna, max_subband_idx] = max(dt_lenna{lev}, [], 3); 
    
    for subband = 1:6
        %for each subband, set to zero all coeffcients that aren't maximal
        %in this subband
 
        temp = dt_lenna{lev}(:,:,subband);
        
        % 1) check we're assigning correctly
        %temp_rm = abs(temp) .* exp(i*angle(temp));
        
        % 2) try random magnitude - as percentage of maximum magnitude
        temp_rm = abs(max_lenna) .* rand(size(temp)) .* exp(i*angle(temp));
        
        temp(max_subband_idx ~= subband) = temp_rm(max_subband_idx ~= subband);
        
        dt_lenna_rm{lev}(:,:,subband) = temp;
        
    end
end

% Now reconstruct the image
lenna_rm = dtwaveifm2(dt_lenna_rm, 'near_sym_b','qshift_b');
figure; imagesc(lenna_rm); axis image; colormap(gray(256));

%%
% Try keeping data from maximal and adjacent sub-bands
dt_lenna_rm2 = dt_lenna;
for lev = 1:4
    
 
    [max_lenna, max_subband_idx] = max(dt_lenna{lev}, [], 3);
    
    for subband = 1:6
        %for each subband, set to zero all coeffcients that aren't maximal
        %in this subband
        
        %for each subband, set to zero all coeffcients that aren't maximal
        %in this subband
        if subband == 6
            sub_up = 1;
        else
            sub_up = subband+1;
        end
        if subband == 1
            sub_down = 6;
        else
            sub_down = subband-1;
        end
        
        temp = dt_lenna{lev}(:,:,subband);
        % 1) check we're assigning correctly
        %temp_rm = abs(temp) .* exp(i*angle(temp));
        
        % 2) try random magnitude - as percentage of maximum magnitude
        temp_rm = abs(max_lenna) .* rand(size(temp)) .* exp(i*angle(temp));
        
        temp_idx = max_subband_idx ~= subband & max_subband_idx ~= sub_up & max_subband_idx ~= sub_down;
        temp(temp_idx) = temp_rm(temp_idx);
        
        dt_lenna_rm2{lev}(:,:,subband) = temp;
        
    end
end

% Now reconstruct the image
lenna_rm2 = dtwaveifm2(dt_lenna_rm2, 'near_sym_b','qshift_b');
figure; imagesc(lenna_rm2); axis image; colormap(gray(256));

%%