clc;
clear all;

pts_per_letter = 300;

% h = text(0,0.5,'ncm','Fontname','elementary sf black','fontsize',250);

%imgpath = 'S:\projects\nailfold\synthesis';
imgpath = 'C:\isbe\nailfold';
pp{1} = ncm_outline(fullfile(imgpath,'b_outline.png'), pts_per_letter);
pp{2} = ncm_outline(fullfile(imgpath,'m_outline.png'), pts_per_letter);
pp{3} = ncm_outline(fullfile(imgpath,'s_outline.png'), pts_per_letter);

pp{2}(:,1) = pp{2}(:,1) + 125;
pp{2}(:,2) = pp{2}(:,2) - 5;
pp{3}(:,1) = pp{3}(:,1) + 220;
pp{3}(:,2) = pp{3}(:,2) + 30;

for i = 1:3
    [widths{i}, dirs{i}, inner{i}, outer{i}] = ...
        generate_edges(pp{i}, 0.005, 42, 1.0);
end

[pp_out, widths_out, inner_out, outer_out, imsz] = ...
    scale_vessel(cat(1,pp{:}), cat(1,widths{:}), ...
                 cat(1,inner{:}), cat(1,outer{:}), 640);

pts = permute(reshape(pp_out', 2, pts_per_letter, 3), [2 1 3]);%unstack(pp_out, [pts_per_letter, 2, 3]);
widths = reshape(widths_out, pts_per_letter, 1, 3);
inner = permute(reshape(inner_out', 2, pts_per_letter, 3), [2 1 3]);
outer = permute(reshape(outer_out', 2, pts_per_letter, 3), [2 1 3]);
dirs = cat(3, dirs{:});

% Trim a few points off the end
n_trim = 0;
pts = pts(1:end-n_trim,:,:);
widths = widths(1:end-n_trim,:,:);
inner = inner(1:end-n_trim,:,:);
outer = outer(1:end-n_trim,:,:);
dirs = dirs(1:end-n_trim,:,:);

max_flow = 4.0;
n_frames = 1000;
samples_per_frame = 1;

for i = 1:3
    [flowmap, mask] = ...
        create_flowmap(pts(:,:,i), widths(:,:,i), dirs(:,:,i), imsz);
    
    flowmap = flowmap * max_flow;

    % Choose cell size and number of cells
    cell_sz = min(widths(:)) * 0.8;
    n_pixels_per_cell = cell_sz^2;
    n_mask_pixels = sum(mask(:)>0);
    n = 10;
    n_cells = ceil(n * n_mask_pixels / n_pixels_per_cell);

    cp{i} = generate_cell_positions(flowmap, mask, ...
                                    pts(:,:,i), widths(:,:,i), ...
                                    n_cells, n_frames+1, ...
                                    [], samples_per_frame);
end
cp = cat(1, cp{:});
    

cloud_add = noiseonf(max(imsz), 1.5);
cloud_add = 32 * normim(cloud_add(1:imsz(1),1:imsz(2)), 'stretch_fixed');

cloud_mult = noiseonf(max(imsz), 1.5);
cloud_mult = 0.75 * normim(cloud_mult(1:imsz(1),1:imsz(2)), 'stretch_fixed');
cloud_mult = 1 - cloud_mult;

outpath = 'C:\isbe\nailfold\synthesis\showcase\ncm_outline\input';
if ~exist(outpath,'dir')
    mkdir(outpath);
else
    delete(fullfile(outpath,'*.png'));
end

f = 1;
for i = 1:n_frames
    img = make_frame(cp(:,:,f:f+samples_per_frame), ...
                     imsz, cell_sz/2, ...
                     [], 140 + cloud_add, 64 * cloud_mult);
	figure(1); clf; hold off; colormap(gray(256));
        imagesc(img); axis image;
    drawnow;

	f = f+samples_per_frame;
    
    filename = sprintf('frame_%04d.png', i);
    imwrite(uint8(img), fullfile(outpath, filename));
end

return

