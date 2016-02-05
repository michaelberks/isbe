% use this function to show the BC coordinate system ontop of the
% breast image
function showCoodSystem(objCS, image)

    if exist('image','var')
        image = processImageMask(image);
        imshow(image,[]); 
    end
    showPoints(objCS)       
    hold on;

    margin = 0.05; % 0.15; %
    for phi = linspace(margin,pi-margin,objCS.resolution)
        s = linspace(0.05,1,objCS.resolution);
        [x, y] = coordSystem2Cart(objCS, s,repmat(phi,1,numel(s))) ;   

        idx = union(find(x<0),find(y<0));
        x(idx)=[]; y(idx)=[];
        plot(y,x,'r-','LineWidth',2);    
    end

    hold off;
end
