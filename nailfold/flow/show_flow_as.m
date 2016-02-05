function img_out = show_flow_as(method, flowmap, varargin)
if (nargin==0 && nargout==0), test(); return; end

switch lower(method)
    case {'uv'},
        % Show as separate horizontal and vertical flowmaps, where green
        % indicates positive (down/right) and red negative (up/left) flow.
        if length(varargin)>0, bg_color = varargin{1}; 
            else bg_color = [0 0 0]; end
        if length(varargin)>1, cmap = varargin{2}; 
            else cmap = redgreen(255); end
        
        % Put U and V side by side
        img = [real(flowmap) nan(size(flowmap,1),1) imag(flowmap)];
            
        % Normalize to range [0..1] where zero is mapped to 0.5
        img = normim(img, 'stretch_fixed');

        % Map to colormap indices, mapping NaN to 1
        img = uint8(1 + ceil(254 * img));

        % Convert to RGB image
        img_out = ind2rgb(img, [bg_color; cmap]);
        
        if (nargout==0)
            image(img_out); axis('image','ij');
        end

    case {'rgb'},
        % Show as RGB map, where hue indicates direction and lightness
        % indicates speed.
        if length(varargin)>0, bgcol = varargin{1}; 
            else bgcol = 1; end
        
        img_out = complex2rgb(flowmap, [],[],[], bgcol);
        
        if (nargout==0)
            image(img_out); axis('image','ij');
        end

    case {'quiver'}
        if length(varargin)>0, nx = varargin{1}; 
            else nx = size(flowmap, 2); end
        if length(varargin)>1, ny = varargin{1}; 
            else ny = size(flowmap, 1); end
            
        if (size(flowmap, 2) < nx) || (nx <= 0)
            nx = size(flowmap, 2); 
        end
        if (size(flowmap, 1) < ny) || (ny <= 0)
            ny = size(flowmap, 1); 
        end
            
        x = round(linspace(1, size(flowmap,2), nx));
        y = round(linspace(1, size(flowmap,1), ny));
        
        if (nargout==0)
            % Show as a grid of flow vectors.
            quiver(x, y, real(flowmap(y,x)), imag(flowmap(y,x)), 2);
        end
        axis image;
        img_out = [];
end
    

%% Test function
function test()
clc;

[xx,yy] = meshgrid(-40:40, -40:40);

rr = sqrt(xx.^2 + yy.^2);

uu = rr .* xx;
vv = rr .* yy;

flowmap = complex(uu,vv);

figure(1); clf;
subplot(3,1,1); show_flow_as('uv', flowmap, [0,0,0], redgreen(255));
subplot(3,1,2); show_flow_as('rgb', flowmap, 1);
subplot(3,1,3); show_flow_as('quiver', flowmap, 20, 20);
