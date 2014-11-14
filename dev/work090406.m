%Code below generates a load of straight vertical lines and measures the
%responses of the DT-CWT - the response should be (equally) maximal in the
%3rd and 4th sub-bands

responses = zeros(128, 8, 6);
x = repmat(linspace(-256, 256, 512), 512, 1);
for halfwidth = 1:128
    
    %pos_line(:, 256-halfwidth:255+halfwidth) = 1;
    sigma2 = (halfwidth^2) / log(2);
    ymax = 1/sqrt(2*pi*sigma2);

    scaling = 1 / ymax;
    pos_line = scaling*exp(-(x.^2 / sigma2)) / sqrt(2*pi*sigma2);
    
    
    dt_pos_line = dtwavexfm2(pos_line, 8);

    for level = 1:8
        r = 256 / 2^level;
        c = 256 / 2^level;

        responses(halfwidth, level, :) = dt_pos_line{level}(r, c, :);

    end
    if ~rem(halfwidth, 32)
        figure; imagesc(pos_line); axis image; colormap gray;
    end
end
%
colors = 'rgbkcmyr';
for band = 1:6
    figure; hold on; title(['Dual-Tree response for band ', num2str(band)]);
    xlabel('Line width');
    ylabel('Magnitude of DT response at centre of line');
    
    for levels = 1:8
        plot(1:128, abs(responses(:,levels,band)), colors(levels));
    end
    
    legend({'Level 1','Level 2','Level 3','Level 4','Level 5','Level 6','Level 7','Level 8'}); 
end
%%
%Now do the same but create lines at different orientations

%We want lines of angle theta, with centre passing through (128,128)
% such a line has equation aX + bY + c = 0 where
% a = sin(theta), b = cos(theta), c = -128*(a + b);
% The distance from any point (x0,y0) to this line is thus:
% |ax0 + by0 + c|

%%
%profile on;

%icp_responses = zeros(128, 180, 8, 3);

x = repmat(1:256, 256, 1);
y = repmat((1:256)', 1, 256);
 
orientations = (1:1:180)*pi/180;
dt_responses = zeros(32, length(orientations), 6, 6);
dt_responses_ri = zeros(32, length(orientations), 6, 6);

%Set central frequencies of bands for interp method
w = [-3 -1; -3 -3; -1 -3; 1 -3; 3 -3; 3 -1]*pi/2.15;

for o = 1:length(orientations)
    a = sin(orientations(o));
    b = cos(orientations(o));
    c = -128*(a + b);
    dx = abs(a*x + b*y + c);
    
    for halfwidth = 1:32
        
        sigma2 = (halfwidth^2) / log(2);
        ymax = 1/sqrt(2*pi*sigma2);

        scaling = 1 / ymax;
        pos_line = scaling*exp(-(dx.^2 / sigma2)) / sqrt(2*pi*sigma2);
        
        dt_pos_line = dtwavexfm2(pos_line, 7);
        dt_pos_line_ri = dtwavexfm2b(pos_line, 7);
        %[ilp_pos_line icp_pos_line] = mb_dual_tree_transform(dt_pos_line);
        
        for level = 1:6
            %r = 128 / 2^level;
            %c = 128 / 2^level;
            
            for band = 1:6
                dt_responses(halfwidth, o, level, band) = ...
                    complex_interp_pixel(dt_pos_line{level}(:, :, band), 128, 128, level, w(band,:));
                dt_responses_ri(halfwidth, o, level, band) = ...
                    complex_interp_pixel(dt_pos_line_ri{level}(:, :, band), 128, 128, level, w(band,:));
            end
            %icp_responses(halfwidth, o, level, :) = icp_pos_line{level}(r, c, :);

        end
        %if ~rem(halfwidth, 32) && ~rem(o, 60)
        %    figure; imagesc(pos_line); axis image; colormap gray;
        %end
    end
end
%profile viewer
% 
colors = 'rgykcmbr';
for halfwidth = 1:32
    
    for level = 1:6
        
        if max(reshape(abs(dt_responses(halfwidth,:,level,:)), 1, [])) > .5
            
            figure;
            subplot(1,2,1); hold on;
            xlabel('Line width');
            ylabel('Magnitude of DT response at centre of line');
            
            for band = 1:6
                plot(orientations, abs(dt_responses(halfwidth,:,level,band)), colors(band));
            end

            legend({'Band 1','Band 2','Band 3','Band 4','Band 5','Band 6'}, 'location', 'southeast');
            title(['Dual-Tree responses for level ', num2str(level), ', halfwidth = ', num2str(halfwidth)]);
            axis([0 pi 0 1]);
            
            subplot(1,2,2); hold on;
            xlabel('Line width');
            ylabel('Magnitude of DT response at centre of line');
            
            for band = 1:6
                plot(orientations, abs(dt_responses_ri(halfwidth,:,level,band)), colors(band));
            end

            legend({'Band 1','Band 2','Band 3','Band 4','Band 5','Band 6'}, 'location', 'southeast');
            title(['Dual-Tree responses for level ', num2str(level), ', halfwidth = ', num2str(halfwidth)]);
            axis([0 pi 0 1]);
        end
        
%         icp_oris = squeeze(icp_responses(halfwidth,:,level,:));
%         icp_oris(imag(icp_oris) < 0) = -icp_oris(imag(icp_oris) < 0);
%         max_icp_oris = max(icp_oris, [], 2);
%         
%         icp_degrees = 180*angle(icp_oris) / pi;
%         max_degrees = 180*angle(max_icp_oris) / pi;
%         subplot(1,2,2); hold on;
%         plot(1:180, 1:180, 'k:');
%         for band = 1:6
%             plot(1:180, icp_degrees(:,band), [colors(band), 'x']); 
%         end
%         plot(1:180, max_degrees, 'LineWidth', 2.0);
%         axis([0,180,0,180]);
%         title('Line orientation vs reported ICP orientation');
    end
end
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
x = repmat(1:256, 256, 1);
y = repmat((1:256)', 1, 256);

halfwidth = 3;
sigma2 = (halfwidth^2) / log(2);
ymax = 1/sqrt(2*pi*sigma2);
scaling = 1 / ymax;

%make lines at angles 23, 157, and vertical reflection of both

%157
orientation = pi*157/180;
a = sin(orientation);
b = cos(orientation);
c = -128*(a + b);
dx = abs(a*x + b*y + c);

line157 = scaling*exp(-(dx.^2 / sigma2)) / sqrt(2*pi*sigma2);
line157_r = flipud(line157);
%23
orientation = pi*23/180;
a = sin(orientation);
b = cos(orientation);
c = -128*(a + b);
dx = abs(a*x + b*y + c);

line23 = scaling*exp(-(dx.^2 / sigma2)) / sqrt(2*pi*sigma2);
line23_r = flipud(line23);

%compute DT
dt_157 = dtwavexfm2(line157, 5);
dt_157_r = dtwavexfm2(line157_r, 5);
dt_23 = dtwavexfm2(line23, 5);
dt_23_r = dtwavexfm2(line23_r, 5);

%Display magnitudes of level 4 coefficients in band 1/6
figure; 
subplot(2,2,1); imagesc(abs(dt_157{4}(:,:,6))); axis image; colormap(gray(256));
subplot(2,2,2); imagesc(abs(dt_157_r{4}(:,:,1))); axis image; colormap(gray(256));
subplot(2,2,3); imagesc(abs(dt_23{4}(:,:,1))); axis image; colormap(gray(256));
subplot(2,2,4); imagesc(abs(dt_23_r{4}(:,:,6))); axis image; colormap(gray(256));

%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Make the diagram explaining this...
figure; hold on; axis equal; axis([120.5 136.5 120.5 136.5]);

[pixels_x pixels_y] = meshgrid(120:137, 120:137);
[l1_x l1_y] = meshgrid(121.5:2:137, 121.5:2:137);
[l2_x l2_y] = meshgrid(122.5:4:137, 122.5:4:137);
[l3_x l3_y] = meshgrid(124.5:8:137, 124.5:8:137);
[l4_x l4_y] = meshgrid(120.5:16:137, 120.5:16:137);

plot(pixels_x(:), pixels_y(:), 'k.');
plot(l1_x(:), l1_y(:), 'ro');
plot(l2_x(:), l2_y(:), 'go');
plot(l3_x(:), l3_y(:), 'bo');
plot(l4_x(:), l4_y(:), 'mo');

legend({'pixels',...
    'level 1 sampling points',...
    'level 2 sampling points',...
    'level 3 sampling points',...
    'level 4 sampling points',...
    }, 'Location', 'EastOutside');
x = [120 137];

%Do stuff for 23 degree line
orientation = pi*23/180;
a = sin(orientation);
b = cos(orientation);
c = -128*(a + b);
y = -(a*x + c) / b;

dxy = abs(a*124.5 + b*124.5 + c);
px = 124.5 + dxy*a;
py = 124.5 + dxy*b;
plot(x, y, 'c');
plot([124.5 px], [124.5 py], 'c:');
text(124.75, 124.75, ['Distance to 23^\o line = ' num2str(dxy)]); 

%Do stuff for 157 degree line
orientation = pi*157/180;
a = sin(orientation);
b = cos(orientation);
c = -128*(a + b);
y = -(a*x + c) / b;

dxy = abs(a*124.5 + b*124.5 + c);
px = 124.5 - dxy*a;
py = 124.5 - dxy*b;
plot(x, y, 'b');
plot([124.5 px], [124.5 py], 'b:');

text(124.75, 124.25, ['Distance to -23^\o line = ' num2str(dxy)]); 

title('Why the discrepancy between DT Coefficients magnitudes in bands 1 and 6?');
