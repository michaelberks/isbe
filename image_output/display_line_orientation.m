colors = hsv(180);

image_in = u_load('C:\isbe\asymmetry_project\data\synthetic_lines\real512\image001.mat');
ori_map = u_load('C:\isbe\asymmetry_project\data\synthetic_lines\real512\results\233902\image001_class.mat');
% s = load('C:\isbe\asymmetry_project\data\synthetic_lines\real512\labels\label001.mat');
% ori_map = s.label_orientation;

%%
f = figure; imagesc(image_in); axis image; colormap(gray(256)); hold on;
ori_map = mod(ori_map-1,180)+1;

while true
    try
        [xi,yi,P] = impixel;
    catch
        break;
    end

    %Plot the coloured line segments
    for ii = 1:length(xi)-1

        deg = ceil(mod(180*atan((yi(ii+1)-yi(ii)) / (xi(ii+1)-xi(ii)))/pi,180));

        plot(xi(ii:ii+1), yi(ii:ii+1), '--', 'color', colors(deg,:));
        
    end
    
    %Plot the orientation estimates as quivers
    [px py ori_estimates] = improfile(ori_map,xi,yi);
    qx = 4*cos(pi*ori_estimates/180);
    qy = 4*-sin(pi*ori_estimates/180);
    for ii = 1:2:length(ori_estimates);
        quiver(px(ii), py(ii), qx(ii), qy(ii), 'color', colors(ceil(ori_estimates(ii)),:));
    end
end
display('finished');
