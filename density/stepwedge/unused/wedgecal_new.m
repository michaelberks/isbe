function wedgevals=wedgecal_new(IMAGE,map,wedgesize)

% find the image size
dimensions = size(IMAGE);

if(wedgesize==1)
    numsteps=25;
elseif(wedgesize==2)
    numsteps=35;
else
    disp('which wedge is this???');
end

disp('select opposite corners of the step wedge area, then press return');
figure(4);
set(4,'Name','Stepwedge Selection');

pointsgiven=1;
while pointsgiven ~= 2

    [swx,swy,P] = impixel(IMAGE,map);
    dimsx=size(swx);
    pointsgiven=dimsx(1);
    
end

% if box goes off edge of image then shift it in a bit:
for ipoint=1:2
    if(swx(ipoint) < 1)
        swx(ipoint)=1;
    elseif(swx(ipoint) > dimensions(2))
        swx(ipoint)=dimensions(2);
    end    
    if(swy(ipoint) < 1)
        swy(ipoint)=1;
    elseif(swy(ipoint) > dimensions(1))
        swy(ipoint)=dimensions(1);
    end
end
    
swimage = IMAGE(swy(1):swy(2),swx(1):swx(2));
dimswimage = size(swimage);
% %cjb find edges and display those too to see if they are informative
% cannyimage=edge(swimage,'canny');
% figure;imshow(cannyimage);

ok=-4;

while (ok<0)

    if(ok<=-3)
        disp('select the edges of two steps, then press return (4 points in all)');
        figure(4);
        set(4,'Name','Stepwedge');
        [swx4,swy4,P] = impixel(swimage,map);

        swx=zeros(2,1);
        swy=zeros(2,1);
    
        swx(1) = (swx4(1)+swx4(2))/2;
        swx(2) = (swx4(3)+swx4(4))/2;
        swy(1) = (swy4(1)+swy4(2))/2;
        swy(2) = (swy4(3)+swy4(4))/2;
    
        line(swx,swy); 
    
        i_min=numsteps+1;
        i_max=0;
        % get step numbers but try to guard against mistakes
        while(i_min<1 | i_min>numsteps | i_max<1 | i_max>numsteps) 
            i_min = input('what number is the first step? ');
            i_max = input('what number is the second step? ');
            %if step numbers accidentally entered wrong way round, don't crash!
            if(i_min>i_max)
                i_temp=i_max;
                i_max=i_min;
                i_min=i_temp;
            end
        end
    
    elseif(ok==-1)
        i_min = i_min+1;
        i_max = i_max+1;
        
        figure(4);
        imshow(swimage,map);
        
    elseif(ok==-2)
        i_min = i_min-1;
        i_max = i_max-1;
        
        figure(4);
        imshow(swimage,map);    
    end
    
    numsteps_selected = i_max - i_min + 1;
    
    % c=improfile;
    % plot(c);

    % draw a line down centre of wedge
    separation=(swy(2)-swy(1))/(numsteps_selected-1);
    offset=(swx(2)-swx(1))/(numsteps_selected-1);
    step_width = sqrt(separation*separation + offset*offset)
    
    wedgevals=zeros(numsteps,1);
    wedgestd=zeros(numsteps,1);
    
    % draw rectangles on each step, representing area to be sampled for grey level
    %cjb put in so debug figures don't overwrite others
    debugfig=20;
    for i=1-i_min:numsteps-i_min
    %   line([swx(1)+i*offset-20 swx(1)+i*offset+20],[swy(1)+i*separation swy(1)+i*separation]);
        figure(4);
        drawrect([swx(1)+i*offset-5 swx(1)+i*offset+5],[swy(1)+i*separation+3 swy(1)+i*separation-3]);
        % cjb mark a couple of steps with step number to aid debug!
        if((i+i_min)==10)
            text(swx(1)+i*offset+5,swy(1)+i*separation,'10');
        elseif((i+i_min)==20)
            text(swx(1)+i*offset+5,swy(1)+i*separation,'20');
        end
        % define indicies of this step's region
        xmin=round(swx(1)+i*offset-5);
        xmax=round(swx(1)+i*offset+5);
        ymin=round(swy(1)+i*separation-3);
        ymax=round(swy(1)+i*separation+3);
        % don't go outside image
        if(ymin<0.5)
            ymin = 1;
        elseif(ymax>dimswimage(1))
            ymax = dimswimage(1);
        end
        if(xmin<0.5)
            xmin = 1;
        elseif(xmax>dimswimage(2))
            xmax = dimswimage(2);
        end
        this_step = swimage(ymin:ymax,xmin:xmax);
        dimensions=size(this_step);
%         %cjb histogram intensities for each step
%         figure(debugfig);
%         imhist(this_step);
%         debugfig=debugfig+1;
        wedgevals(i+i_min) = mean(double(this_step(:)));
        wedgestd(i+i_min) = std(double(this_step(:)));
        if(wedgestd(i+i_min) > 7.5)
            figure(debugfig);
            imhist(this_step);
            debugfig=debugfig+1;
            disp('step has large std - look at intensity histogram');
            i+i_min
            disp('press enter to continue...');
            pause;
        end
        %%%%%%% show histogram of grey levels on this step %%%%%%%%%%%%%%%
        %figure(5);
        %hist(double(this_step(:)));
        %figure(4);
        %disp(mean(double(this_step(:))))
        %disp(std(double(this_step(:))))
        %pause;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %mask=roipoly(swimage,[swx(1)+i*offset-20 swx(1)+i*offset+20 swx(1)+i*offset+20 swx(1)+i*offset-20],[swy(1)+i*separation+10 swy(1)+i*separation+10 swy(1)+i*separation-10 swy(1)+i*separation-10]);
        %pixval(i)=sum(double(swimage).*double(mask));
    end
    
    disp('Are these regions ok?');
    disp('enter "-3" for no, define points again');
    disp('enter "-2" for no, nudge up a step');
    disp('enter "-1" for no, nudge down a step');
    disp('enter "0" for yes');
    disp('enter "2" for fix needed');
    ok=input(':');

end

while(ok==2)
    
    figure(8);
    plot(wedgevals);
    
    bad_first = input('what is first bad step? ');
    bad_last = input('what is last bad step? ');
    
    if(bad_first == 1)
        wedgevals(bad_first:bad_last) = wedgevals(bad_last+1);
        wedgestd(bad_first:bad_last) = wedgestd(bad_last+1);
    elseif(bad_last == numsteps)
        wedgevals(bad_first:bad_last) = wedgevals(bad_first-1);
        wedgestd(bad_first:bad_last) = wedgestd(bad_first-1);
    else
        temp1 = [1:bad_last-bad_first+1];
        temp2 = temp1.*(wedgevals(bad_last+1)-wedgevals(bad_first-1))/(bad_last-bad_first+2) + wedgevals(bad_first-1);
        wedgevals(bad_first:bad_last)=temp2;    
    end
    figure(8);
    plot(wedgevals);
    
    % check if this is OK
    disp('Is this OK?');
    disp('enter "1" for yes');
    disp('enter "2" for no, try again');
    ok = input(':');
   
end

% first check there isn't a problem with the first step -
% if this goes too close to edge of film, then it may not be as dark as it should

if(wedgevals(1) < wedgevals(2))
    disp('problem with first step - not as dark as it should be');
    wedgevals(1) = wedgevals(2);
    wedgestd(1) = wedgestd(2);
end

figure(debugfig);
plot(wedgevals,'bo');
xlabel('step number');
ylabel('mean intensity');
title('stepwedge intensity vs step number');
debugfig=debugfig+1;
% correct steps which have unphysical pixel value (eg. due to scatter/edge effects or stepwedge obstruction)
% i.e. as thickness of wedge increases, pixel value can only decrease.
% cjb commented this out because it does really really stupid things.
% 28/6/04
halfwaystep=round((numsteps-1)/2);
for i=halfwaystep:numsteps-1
    if(wedgevals(i)<wedgevals(i+1))
        disp('step above this has intensity higher - setting same');
        i
        wedgevals(i+1)=wedgevals(i);
        wedgestd(i+1) = wedgestd(i);
    end
end
for i=halfwaystep:-1:2
    if(wedgevals(i-1)<wedgevals(i))
        disp('step below this has intensity lower - setting same');
        i
        wedgevals(i-1)=wedgevals(i);
        wedgestd(i-1) = wedgestd(i);
    end
end

% for i=1:numsteps-1
%     if(wedgevals(i)<wedgevals(i+1))
%         disp('step above this has intensity higher - setting same');
%         i
%         wedgevals(i+1)=wedgevals(i);
%         wedgestd(i+1) = wedgestd(i);
%     end
% end
hold on;
plot(wedgevals,'r+');
legend('before correction','after correction');
hold off;

for i=1:numsteps
    if(wedgestd(i) > 7.5)
        disp('possible large std dodgyness on step:');
        disp(i);
        disp('press enter to continue...');
        pause
    end
end