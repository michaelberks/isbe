function breast_thickness = edge_profile(x_b,mask_b,edgex,edgey,leftright,resize_factor)
%leftright=0:right, 1: left

dims=size(mask_b);

if(leftright==0)
    mask_b=flipdim(mask_b,2);
    % flip edge points too
    edgey=edgey*-1. + dims(2) + 1;
    %flip x_b cos mask has been flipped
    x_b=flipdim(x_b,2);
end

masksize=size(mask_b);
xbsize=size(x_b);
figure(51);imshow(mask_b,[]);

figure(52);imshow(x_b,[]);
figure(53);plot(edgex,edgey);

% now have edge, mask and thickness images all same size.
breast_thickness=x_b;
breast_thickness(x_b>0)=0;

% get points +/-10 of the centre point
localedge=zeros(2,21);
%loop over edge points, leaving out 10
numedgepoints=size(edgex);
minx=min(edgex);
maxx=max(edgex);
%for ix=10:3:numedgepoints-10 % try modifying these parameters
for ix=10:3:numedgepoints-10 % try modifying these parameters
    if(edgex(ix)<minx+20 || edgex(ix)>maxx-20 || edgey(ix)<20 || edgey(ix)>dims(2)-20)
%         disp('Too close to edge of image');
    else
        centrex=edgex(ix);
        centrey=edgey(ix);
        pointnum=1;
        for jx=ix-10:ix+10
            localedge(1,pointnum)=edgex(jx);
            localedge(2,pointnum)=edgey(jx);
            pointnum=pointnum+1;
        end
        [p,s]=polyfit(localedge(2,:),localedge(1,:),1);
        figure(518);
        h=line(localedge(2,:),localedge(1,:));
        set(h,'color','red');
        
        % gradient is p(1) so normal gradient is -1/p(1)
        % want line passing through centrex but perp to line - 
        % yspatial=(-1/p(1))xspatial + (centrex + edge_points(centrex)/p(1)
        % THIS BIT IS OK
        normgrad=-1/p(1);
        
        xnorm=100;
        ynorm=normgrad*xnorm+(centrex-centrey*normgrad);
        h3=line([centrey xnorm],[centrex ynorm]);
        set(h3,'color','blue');
        
        % get thickness at the centre point, divide by 2 to get radius - need to
        % correct to pixels not mm!!!
        centrethick=double(x_b(centrex,centrey))*1e-3*resize_factor/(50e-6);
        centrethick=centrethick/2.;
        % calculate point on the normal which is centrethick away from the edge
        delx=centrethick*cos(atan(normgrad));
        dely=centrethick*sin(atan(normgrad));
        endpointx=double(centrey)-delx;
        endpointy=centrex-dely;
        
        h1=line([centrey endpointx],[centrex endpointy]);
        set(h1,'color','green');
        
        length=sqrt(delx^2+dely^2);
        
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
        numprofpoints=round(length)+1;
        steplength=length/numprofpoints;
        for j=1:numprofpoints
            %find point position
        delx=j*steplength*cos(atan(normgrad));
        dely=j*steplength*sin(atan(normgrad));
        endpointx=round(double(centrey)-delx);
        endpointy=round(double(centrex)-dely);
            if(endpointx>0 && endpointy>0)
                %get distance between profile point and edge in pixels
                edgedist=delx*delx+dely*dely;
                edgedist=sqrt(edgedist);
                %calculate thickness estimate based on semicircular profile in pixels
                breastthick=double(x_b(endpointy,endpointx))*1e-3*resize_factor/(50e-6);
                edgethick=2.*sqrt(abs((breastthick-edgedist)*edgedist));
                %compare with thickness of breast - if <, enter into profile, if >=, stop
                if(edgethick<breastthick)
                    breast_thickness(endpointy,endpointx)=edgethick;
                end
            end
        end
    end
end

%figure(56);imshow(breast_thickness,[]);

% now fill in gaps
gapmask=breast_thickness;
gapmask(:,:)=0;
gapmask(breast_thickness==0)=1;
%figure(57);imshow(gapmask,[]);
edgesetmask=breast_thickness;
edgesetmask(:,:)=0;
edgesetmask(breast_thickness>0)=1;
%figure(58);imshow(edgesetmask,[]);
avefilter=ones(5);
% find number of set pixels in neighbourhood 
neighbours=imfilter(edgesetmask,avefilter);
%set zero neighbours to 1 to avoid divide by zero error - don't use these
%pixels anyway
neighbours(neighbours==0)=1;
%figure(511);imshow(neighbours,[]);
%fill in gaps by using averaging filter
gapthick=imfilter(double(breast_thickness),avefilter);
%correct for number of neighbours used to calculate the average
gapthick=double(gapthick)./double(neighbours);
gapthick=gapthick.*double(gapmask);
%figure(59);imshow(gapthick,[]);
breast_thickness=double(breast_thickness)+double(gapthick);
%figure(510);imshow(breast_thickness,[]);
colormap hot;
breast_thickness_smoo=medfilt2(breast_thickness);
%figure(512);imshow(breast_thickness_smoo,[]);
colormap hot;
%now add rest of thickness map - remember to convert edge thickness back to
%mm for consistency
breast_thickness_smoo=double(breast_thickness_smoo)*(50e-6)/(1e-3*resize_factor);
% remove pixels outside the breast mask
overlapim=breast_thickness_smoo;
overlapim(mask_b==0)=0;
figure(13);imshow(overlapim,[]);

%flip mask and edge thickness back
if(leftright==0)
    mask_b=flipdim(mask_b,2);
    overlapim=flipdim(overlapim,2);  
    x_b=flipdim(x_b,2);
end


% combine with plate separation (x_b) image
sepim=double(x_b).*double(mask_b)/max(max(double(mask_b)));
sepim(overlapim>0)=0;
%figure(514);imshow(sepim,[]);
all_breast_thickness=sepim+overlapim;
figure(515);imshow(all_breast_thickness,[]);
imwrite(all_breast_thickness,'map.bmp','bmp');

breast_thickness=all_breast_thickness;
disp('max breast_thickness in edge_profile');
max(max(breast_thickness))
disp('min breast_thickness in edge_profile');
min(min(breast_thickness))
