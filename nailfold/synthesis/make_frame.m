function [img, flow_img] = make_frame(cell_positions, imsz, ...
                                      blob_hw, mask, ...
                                      background, contrast, ...
                                      n_exposed)

if (nargin==0 && nargout==0), test_script(); return; end
                
if ~exist('mask','var'), mask = ones(imsz); end
if ~exist('background','var'), background = 120; end
if ~exist('contrast','var'), contrast = 10; end
if ~exist('n_exposed','var'), n_exposed = inf; end

n_exposed = max(min(n_exposed, size(cell_positions,3)-1), 1);

% Define the blob
x = linspace(-3, 3, 2*blob_hw+1);
dx = x(2)-x(1);

n_cells = size(cell_positions,1);

%% Compute observed image
img = ones(imsz);

% Update synthetic image
for i = 1:n_exposed
    for c = 1:n_cells
        p = cell_positions(c, 1:2, i);
        [rrng,crng] = blob_pixels(p, blob_hw, size(img));
        
        dr = (rrng - p(2)) * dx;
        dc = (crng - p(1)) * dx;
        gr = exp(-0.5 * dr.*dr);
        gc = exp(-0.5 * dc.*dc);
        blob = (gr - min(gr(:)))' * (gc - min(gc(:)));

        opacity = 0.25 * cell_positions(c, 3, i);
%         opacity = 1;

        % Darken the image by the value opacity
        blob = (1 - opacity*blob/n_exposed);
        img(rrng, crng) = img(rrng, crng) .* blob;
    end
end

% % Mask to tidy up edges and downweight thin vessels
% img = 1-((1-img).*mask);
% 
% % Smooth the image a little.
% persistent g;
% if isempty(g)
%     s = 4;
%     x = linspace(-s*4,s*4, s*8);
%     g = exp(-x.*x/(2*s*s));
%     g = g / sum(abs(g));
% end
%
% % Convolve with 1-img to avoid boundary effects
% img = conv2(g, g, 1-img, 'same');

img = 1-img;

img = ... % background value/image
      background - ... 
      ... % scaled vessel image
      contrast .* img; 
  
use_speckle_noise = true;
if use_speckle_noise
    g_min = min(img(:));
    g_max = max(img(:));
    img2 = (img - g_min) / (g_max-g_min);
    img2 = imnoise(img2, 'speckle', 0.01);
    img = img2*(g_max-g_min) + g_min;
end

if (nargout==1), return; end
  
%% Compute observed flow
flow_img = complex(zeros([imsz,2]), zeros([imsz,2]));

if (size(cell_positions,3) > 1)
    displacements = cell_positions(:,:,end) - cell_positions(:,:,1);

    for c = 1:n_cells
        reborn_cell = (displacements(c, 4) ~= 0);

        if (nargout > 1) && ...
           (size(cell_positions, 3) > 1) && ...
           (~reborn_cell)

            p = cell_positions(c, 1:2, i);
            [rrng,crng] = blob_pixels(p, blob_hw, size(img));
            
            % Update the accumulated flow vector
            flowvec = complex(displacements(c, 1), displacements(c, 2));
            flow_img(rrng, crng, 1) = flow_img(rrng, crng, 1) + flowvec;

            % Increase the count for these pixels
            flow_img(rrng, crng, 2) = flow_img(rrng, crng, 2) + 1;
        end
    end
end


%% Determine which rows and columns are occupied by the blob at p
function [rrng, crng] = blob_pixels(p, blob_hw, imsz)

rrng = floor(p(2)-blob_hw):ceil(p(2)+blob_hw);
rvalid = ((rrng > 0) & (rrng < imsz(1)));
rrng = rrng(rvalid);

crng = floor(p(1)-blob_hw):ceil(p(1)+blob_hw);
cvalid = ((crng > 0) & (crng < imsz(2)));
crng = crng(cvalid);


%% Test script
function test_script()
clc;

s = randn('state');
disp(uint32(s));

f_profile = false;
if f_profile, profile clear; profile on; end

pts = generate_path(500);

[widths, dirs, inner, outer] = ...
    generate_edges(pts, [], 10, 1.6);

[pts, widths, inner, outer, imsz] = ...
    scale_vessel(pts, widths, inner, outer, 160);

[flowmap, mask] = ...
    create_flowmap(pts, widths, dirs, imsz);

max_flow = 40.0;
flowmap = flowmap * max_flow;

n_cells = 20;
n_frames = 200;
frame_skip = 8;
cell_positions = generate_cell_positions(flowmap, mask, ...
                                         pts, widths, ...
                                         n_cells, n_frames+1, ...
                                         [], [], ...
                                         frame_skip);
figure(1); clf;
subplot(1,2,1);
    animate_cell_positions(cell_positions(:,:,1:frame_skip:end), imsz);
                                       
blob_hw = 5; %min(widths);
img = make_frame(cell_positions(:,:,1:1+frame_skip), imsz, blob_hw, ...
                 [], 200, 64);

if f_profile, profile off; profile report; end

f = 1;
for i = 1:n_frames
    img = make_frame(cell_positions(:,:,f:f+frame_skip), imsz, blob_hw, ...
                     [], 200, 64);
	figure(1); colormap(gray(256));
    subplot(1,2,2); cla; 
        image(img); axis('image');
    drawnow;

	f = f+frame_skip;
end
    
return