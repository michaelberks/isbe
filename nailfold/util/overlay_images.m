imgpath = 'U:\projects\nailfold\synthesis\showcase\zoom\ncm_video';

overlay = double( imread(fullfile(imgpath, 'rgb_overlay.png')) );
legend = double( imread(fullfile(imgpath, 'rgb_legend.png')) );
legend_mask = repmat((mean(legend, 3) ~= 255), [1,1,3]);

d = dir(fullfile(imgpath, 'frame_*.png'));

outpath = fullfile(imgpath, 'overlaid');
if ~exist(outpath,'dir')
    mkdir(outpath);
else
%     delete(fullfile(outpath,'frame_*.png'));
end

alpha = linspace(1,0,ceil(length(d)/2));
alpha = 0.5 * [alpha 1-alpha];

for i = 1:length(d)
    img = double( imread(fullfile(imgpath, d(i).name)) );
    
    rgb = (1-alpha(i)) * img + ...
             alpha(i)  * overlay;
      
    rgb(1:size(legend,1), 1:size(legend,2), :) = ...
        legend_mask .* legend + ...
        ~legend_mask .* rgb(1:size(legend,1), 1:size(legend,2), :);
   
    filename = sprintf('frame_%04d.png', i);
    imwrite(uint8(rgb), fullfile(outpath, filename));
    
%     figure(1); clf;
%         image(uint8(rgb));
%         axis('image');
%     drawnow;
end