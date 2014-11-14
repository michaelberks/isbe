function [ image, anon_map ] = anon_white( image )
%ANON_WHITE : Turn anon. region white instead of black to make it good for cross correlation and map the area
% 
%
% Inputs:
%			image				image to search in [2D array]
%
% Outputs:
%			image				new image with changed anon region [2D array]
%			anon_map 			map of the anon region (0 over anon region, 1 everywhere else) [2D array]
%
% Example:
%
% Notes:	
%			
% See also:
%
% Created: 18-06-2013
% Author: Euan Allen
% Email : euan.allen@gmail.com

N = 1000;%how many pixels into image to check
[Y X] = size(image);
anon_map = ones(size(image)); %to track when anon is
count = 0;
limit = 0.25; %fraction of image edge that has to be continuously black to consider it the anon region

%left side
for i = 1:N %search N edges in
    temp = image(:,i); %take edge strip as temp
    for ii = 1:Y %search along strip
        if temp(ii) ==0&&ii~=Y %if its a black pixel count or end
            count = count+1;
        elseif count~=0 %if its not a black pixel and we've been counting
            if count>limit*Y %if consecutive black region is large enough
                %image(ii-count:ii-1,i)=0; %keep black
                image(ii-count:ii-1,i)=2^16; %make white
                anon_map(ii-count:ii-1,i)=0; %make map of the region 0
            end
            count = 0; %reset
        end
    end
end
clear temp

%right side
for i = X-N-1:X 
    temp = image(:,i);
    for ii = 1:Y
        if temp(ii) ==0&&ii~=Y 
            count = count+1;
        elseif count~=0 
            if count>limit*Y 
                image(ii-count:ii-1,i)=2^16;
                anon_map(ii-count:ii-1,i)=0;
            end
            count = 0; %reset
        end
    end
end
clear temp

%top side
for i = 1:N 
    temp = image(i,:);
    for ii = 1:X
        if temp(ii) ==0&&ii~=X 
            count = count+1;
        elseif count~=0 
            if count>limit*X 
                image(i,ii-count:ii-1)=2^16;
                anon_map(i,ii-count:ii-1)=0;
            end
            count = 0; %reset
        end
    end
end
clear temp

%bottom side
for i = Y-N-1:Y
    temp = image(i,:);
    for ii = 1:X
        if temp(ii) ==0&&ii~=X 
            count = count+1;
        elseif count~=0 
            if count>limit*X 
                image(i,ii-count:ii-1)=2^16;
                anon_map(i,ii-count:ii-1)=0;
            end
            count = 0; %reset
        end
    end
end
end

