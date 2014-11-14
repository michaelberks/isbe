function animate_cell_positions(cell_positions, imsz, f_show_vector)
if (nargin==0 && nargout==0), test(); return; end

if ~exist('imsz','var') || isempty(imsz), imsz = max(max(cell_positions(:,[2,1],:),[],3),[],1) * 1.1; end
if ~exist('f_show_vector','var') || isempty(f_show_vector), f_show_vector = false; end

nCells = size(cell_positions, 1);
nFrames = size(cell_positions, 3);

for i = 1:nFrames
    cla; hold on;
    
    if f_show_vector && (i > 2)
        cell_x = cell_positions(:, 1, i-1:i);
        cell_y = cell_positions(:, 2, i-1:i);
        plot(cell_x(:,:)', cell_y(:,:)', 'r-');
    end
    plot(cell_positions(:,1,i), cell_positions(:,2,i), 'k.', ...
         'markersize', 16);
    axis('equal','ij', [1,imsz(2),1,imsz(1)]);
    set(gca,'box','on','xticklabel',[],'yticklabel',[]);
    drawnow;
end


%% Test script    
function test()
clc;

nCells = 100;
nFrames = 100;

cell_positions = zeros(nCells, 2, nFrames);
cell_positions(:,:,1) = rand(nCells, 2);
for i = 2:nFrames
    cell_positions(:,:,i) = cell_positions(:,:,i-1) + ...
                            0.01*randn(nFrames, 2);
end

cell_positions = 160*cell_positions;
animate_cell_positions(cell_positions, 160*[1,1]);
