function breast_thickness = edge_profile(x_b, mask_b, edgex, edgey, leftright, resize_factor, debug_mode)
%leftright=0:right, 1: left

if nargin < 7
    debug_mode = 0;
end

dims=size(mask_b);

if(leftright==0)
    mask_b=flipdim(mask_b,2);
    % flip edge points too
    edgex=edgex*-1. + dims(2) + 1;
    %flip x_b cos mask has been flipped
    x_b=flipdim(x_b,2);
end


% now have edge, mask and thickness images all same size.
breast_thickness = x_b;
breast_thickness(x_b>0) = 0;

% get points +-10 of the centre point
localedge = zeros(2,21);
%loop over edge points, leaving out 10
numedgepoints = length(edgey);
miny = min(edgey);
maxy = max(edgey);

if debug_mode
    f1 = figure; hold on;
end
for ix = 10:1:numedgepoints-10
    
    if(edgey(ix)<miny+20 || edgey(ix)>maxy-20 || edgex(ix)<20 || edgex(ix)>dims(2)-20)
        %disp('Too close to edge of image');
    else
        centrex = edgey(ix);
        centrey = edgex(ix);
       
        [p] = polyfit(edgex(ix-10:ix+10), edgey(ix-10:ix+10),1);

        % gradient is p(1) so normal gradient is -1/p(1)
        % want line passing through centrex but perp to line - 
        % yspatial=(-1/p(1))xspatial + (centrex + edge_points(centrex)/p(1)
        % THIS BIT IS OK
        normgrad = -1/p(1);
        
        xnorm = 100;
        ynorm = normgrad*xnorm+(centrex-centrey*normgrad);
        
        
        
        % get thickness at the centre point, divide by 2 to get radius - need to
        % correct to pixels not mm!!!
        centrethick = double(x_b(centrex,centrey))*1e-3*resize_factor/(44e-6);
        centrethick = centrethick/2.;
        % calculate point on the normal which is centrethick away from the edge
        delx = centrethick*cos(atan(normgrad));
        dely = centrethick*sin(atan(normgrad));
        endpointx = double(centrey)-delx;
        endpointy = centrex-dely;
               
        norm_length = sqrt(delx^2+dely^2);
        
        % Now I have a decent line, get image profile to fill with thicknesses
        % do on image read in, remembering coordinates profx profy are SPATIAL
        % figure(2);
        %[profx,profy,profile]=improfile(x_b,double([centrey endpointx]),double([centrex endpointy]));
        % figure(5);
        % plot(profile);
        %numprofpoints=size(profx);
        
        
%         figure(518);
%         h2=line(profx,profy);
%         set(h2,'color','blue');
%         set(h2,'linestyle',':');

%        for j=1:numprofpoints(1)     

        % calculate profile points in loop instead.  Length is in pixels so
        % use that as an estimate of how many points to use
        numprofpoints = round(norm_length)+1;
        steplength = norm_length/numprofpoints;
        for j = 1:numprofpoints
            %find point position
            delx = j*steplength*cos(atan(normgrad));
            dely = j*steplength*sin(atan(normgrad));
            endpointx = round(double(centrey)-delx);
            endpointy = round(double(centrex)-dely);
            if (endpointx>0 && endpointy>0 && endpointy <= dims(1))
                %get distance between profile point and edge in pixels
                edgedist = delx*delx+dely*dely;
                edgedist = sqrt(edgedist);
                %calculate thickness estimate based on semicircular profile in pixels
                breastthick = double(x_b(endpointy,endpointx))*1e-3*resize_factor/(44e-6);
                edgethick = 2.*sqrt(abs((breastthick-edgedist)*edgedist));
                %compare with thickness of breast - if <, enter into profile, if >=, stop
                if (edgethick < breastthick) %#ok
                    breast_thickness(endpointy,endpointx)=edgethick;
                end
            end
        end
        
        if debug_mode
            figure(f1);
            plot(localedge(2,:),localedge(1,:), 'r');
            plot([centrey xnorm],[centrex ynorm], 'b');
            plot([centrey endpointx],[centrex endpointy], 'g');
        end
        
    end
end


% now fill in gaps
gapmask=breast_thickness;
gapmask(:,:)=0;
gapmask(breast_thickness==0) = 1;

edgesetmask = breast_thickness;
edgesetmask(:,:) = 0;
edgesetmask(breast_thickness>0) = 1;

avefilter = ones(5);

% find number of set pixels in neighbourhood 
neighbours=imfilter(edgesetmask,avefilter);

%set zero neighbours to 1 to avoid divide by zero error - don't use these
%pixels anyway
neighbours(neighbours==0)=1;

%fill in gaps by using averaging filter
gapthick = imfilter(double(breast_thickness),avefilter);

%correct for number of neighbours used to calculate the average
gapthick = double(gapthick)./double(neighbours);
gapthick = gapthick.*double(gapmask);

breast_thickness = double(breast_thickness)+double(gapthick);
breast_thickness_smoo = medfilt2(breast_thickness);

%now add rest of thickness map - remember to convert edge thickness back to
%mm for consistency
breast_thickness_smoo = double(breast_thickness_smoo)*(44e-6)/(1e-3*resize_factor);

% remove pixels outside the breast mask
overlapim = breast_thickness_smoo;
overlapim(mask_b==0) = 0;

%flip mask and edge thickness back
if(leftright==0)
    mask_b = flipdim(mask_b,2);
    overlapim = flipdim(overlapim,2);  
    x_b = flipdim(x_b,2);
end


% combine with plate separation (x_b) image
sepim = double(x_b).*double(mask_b)/max(max(double(mask_b)));
sepim(overlapim>0) = 0;

breast_thickness = sepim+overlapim;

%MB only show in debug_mode
if debug_mode
    figure; imshow(breast_thickness,[]);
    
    disp(['max breast_thickness in edge_profile = ', num2str(max(breast_thickness(:)))]);
    disp(['min breast_thickness in edge_profile = ', num2str(min(breast_thickness(:)))]);

end

%MB I don't think we want to do this
%imwrite(breast_thickness, 'map.bmp', 'bmp');





