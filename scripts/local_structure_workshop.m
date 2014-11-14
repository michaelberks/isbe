%--------------------------------------------------------------------------
% Some useful Matlab code for generating test images, and applying methods
% to detect the local structure in these images
%--------------------------------------------------------------------------

%%
%--------------------------------------------------------------------------
%1: Generate Gaussian bars of varying widths and orientations in 256x256 images
%   - save the DT-CWT response at the centre of the line
%   - can also save ILP/ICP responses
%   - could add in other loval structure methods (e.g. monogenic signal)

%Allocate space for Dt/ICP responses
dt_responses = zeros(128, 180, 8, 6);
icp_responses = zeros(128, 180, 8, 6);

%Generate x-y coordinate matrices
x = repmat(1:256, 256, 1);
y = repmat((1:256)', 1, 256);

%Choose the orientations of the lines (1-180 degrees here)
orientations = (1:180)*pi/180;

%for each orientation
for o = 1:length(orientations)
    
    %work out the distance of each pixel 
    a = sin(orientations(o));
    b = cos(orientations(o));
    c = -128*(a + b);
    dx = abs(a*x + b*y + c);
    
    %for each width of line (defined as the width of the Gaussian bar at
    %half max height
    for halfwidth = 1:128
        
        %Compute corresponding Gaussian sigma and max_height
        sigma2 = (halfwidth^2) / log(2);
        ymax = 1/sqrt(2*pi*sigma2);
        
        %Compute scaling to have max height == 1
        scaling = 1 / ymax;
        
        %Generate Gaussian bar
        pos_line = scaling*exp(-(dx.^2 / sigma2)) / sqrt(2*pi*sigma2);
        
        %Compute local structure responses
        dt_pos_line = dtwavexfm2(pos_line, 7);
        [ilp_pos_line icp_pos_line] = mb_dual_tree_transform(dt_pos_line);
        
        %For each level of decomposition, save response at bar centre 
        for level = 1:6
            r = 128 / 2^level;
            c = 128 / 2^level;

            dt_responses(halfwidth, o, level, :) = dt_pos_line{level}(r, c, :);
            icp_responses(halfwidth, o, level, :) = icp_pos_line{level}(r, c, :);

        end
        %if ~rem(halfwidth, 32) && ~rem(o, 60)
        %    figure; imagesc(pos_line); axis image; colormap gray;
        %end
    end
end

%%
% 2. Generate bars as in 1) but with pairs of bars crossing - will just do
% this for a fixed halfwidth
%   - we'll just display the results for now but could save as before...

%Generate x-y coordinate matrices
x = repmat(1:256, 256, 1);
y = repmat((1:256)', 1, 256);

%Choose halfwidth
halfwidth = 5;

%Compute corresponding Gaussian sigma and max_height
sigma2 = (halfwidth^2) / log(2);
ymax = 1/sqrt(2*pi*sigma2);

%Compute scaling to have max height == 1
scaling = 1 / ymax;

%Choose the orientations of the lines (18:18:180 degrees here)
orientations = (18:18:180)*pi/180;

%Create horizontal bar - a=0, b=1, c=-128, dx = abs(y-120.5)
dx = abs(y - 120.5);
bar_horizontal = scaling*exp(-(dx.^2 / sigma2)) / sqrt(2*pi*sigma2);

colors = 'rgymcb';

%for each orientation
for o = 1:length(orientations)
    
    %work out the distance of each pixel 
    a = sin(orientations(o));
    b = cos(orientations(o));
    c = -120.5*(a + b);
    dx = abs(a*x + b*y + c);
    
    %generate oriented bar
    bar_oriented = scaling*exp(-(dx.^2 / sigma2)) / sqrt(2*pi*sigma2);
    
    %make bar cross
    bar_cross = max(bar_oriented, bar_horizontal);

    %Compute local structure responses
    dt = dtwavexfm2(bar_cross, 7);
    [ilp icp] = mb_dual_tree_transform(dt);
    
    %Plot the ILP/ICP responses
    
    for level = 4

        st = .5 - 2^(level - 1);
        pts = st + (2^level)*(1:2^(8-level));
        
        sf_ilp = (2^level) / max(abs(ilp{level}(:)));
        sf_icp = (2^level) / max(abs(icp{level}(:)));
        
        f_ilp = figure('Name', ['ILP level ', num2str(level)]);
        imagesc(bar_cross); axis image; colormap(gray(256)); hold on;        
        for band = 1:6
            quiver(pts, pts, sf_ilp*real(ilp{level}(:,:,band)), -sf_ilp*imag(ilp{level}(:,:,band)), 0, colors(band));
        end
        
        f_icp = figure('Name', ['ICP level ', num2str(level)]);
        imagesc(bar_cross); axis image; colormap(gray(256)); hold on;        
        for band = 1:6
            quiver(pts, pts, sf_icp*real(icp{level}(:,:,band)), -sf_icp*imag(icp{level}(:,:,band)), 0, colors(band));
        end
 
    end
end
%%
%--------------------------------------------------------------------------
% 3. Build a Gaussian bar/edge and display the interpolated dual-tree and ILP
% coefficients at each pixel

x = repmat(1:256, 256, 1);
y = repmat((1:256)', 1, 256);

%Choose halfwidth
halfwidth = 5;

%Compute corresponding Gaussian sigma and max_height
sigma2 = (halfwidth^2) / log(2);
ymax = 1/sqrt(2*pi*sigma2);

%Compute scaling to have max height == 1
scaling = 1 / ymax;

%work out the distance of each pixel 
a = sin(-pi/6);
b = cos(-pi/6);
c = -128*(a + b);
dx = abs(a*x + b*y + c);
    
%generate oriented bar
gaussian_bar = scaling*exp(-(dx.^2 / sigma2)) / sqrt(2*pi*sigma2);
display_full_tree_info(gaussian_bar);
%phase_map = dt_phase_congruency(gaussian_bar, 2, 5);
%figure; imagesc(complex2rgb(max(phase_map, [], 3))); axis image; colormap(gray(256));

gaussian_bar2 = 1 - scaling*exp(-(dx.^2 / sigma2)) / sqrt(2*pi*sigma2);
display_full_tree_info(gaussian_bar2);
%phase_map = dt_phase_congruency(gaussian_bar2, 2, 5);
%figure; imagesc(complex2rgb(max(phase_map, [], 3))); axis image; colormap(gray(256));
%
dx2 = (a*x + b*y + c);
dx2(dx2 < 0) = 0;
gaussian_edge = scaling*exp(-(dx2.^2 / sigma2)) / sqrt(2*pi*sigma2);
display_full_tree_info(gaussian_edge);
%phase_map = dt_phase_congruency(gaussian_edge, 2, 5);
%figure; imagesc(complex2rgb(max(phase_map, [], 3))); axis image; colormap(gray(256));
%
dx3 = (a*x + b*y + c);
dx3(dx3 > 0) = 0;
gaussian_edge2 = scaling*exp(-(dx3.^2 / sigma2)) / sqrt(2*pi*sigma2);

display_full_tree_info(gaussian_edge2);
%phase_map = dt_phase_congruency(gaussian_edge2, 2, 5);
%figure; imagesc(complex2rgb(max(phase_map, [], 3))); axis image; colormap(gray(256));

%%
%--------------------------------------------------------------------------
%5: Generate Gaussian bars of varying widths and orientations in 256x256 images
%   - save the DT-CWT response at the centre of the line
%   - can also save ILP/ICP responses
%   - could add in other loval structure methods (e.g. monogenic signal)

%Allocate space for Dt/ICP responses
dt_responses = zeros(64, 90, 6, 7);


%Generate x-y coordinate matrices
x = repmat(1:256, 256, 1);
y = repmat((1:256)', 1, 256);

%Choose the orientations of the lines (1-180 degrees here)
orientations = (2:2:180)*pi/180;

%for each orientation
for o = 1:length(orientations)
    
    %work out the distance of each pixel 
    a = sin(orientations(o));
    b = cos(orientations(o));
    c = -128*(a + b);
    dx = abs(a*x + b*y + c);
    
    %for each width of line (defined as the width of the Gaussian bar at
    %half max height
    for halfwidth = 1:64
        
        %Compute corresponding Gaussian sigma and max_height
        sigma2 = (halfwidth^2) / log(2);
        ymax = 1/sqrt(2*pi*sigma2);
        
        %Compute scaling to have max height == 1
        scaling = 1 / ymax;
        
        %Generate Gaussian bar
        pos_line = scaling*exp(-(dx.^2 / sigma2)) / sqrt(2*pi*sigma2);
        
        %Compute local structure responses
        dt_pos_line = dtwavexfm2(pos_line, 7);
        
        %Compute full tree
        ft_pos_line = dt_to_full_image(dt_pos_line);
        clear dt_pos_line;
        
        
        %For each subband of decomposition, save response at bar centre
        dt_responses(halfwidth, o, :, :) = ft_pos_line(128,128,:,:);
    end
end
save C:\isbe\dev\background\orientations\dt_responses_int.mat dt_responses
