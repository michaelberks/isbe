%function [phi_x_out] = Phasefield2D_binary(num)
num = 33;

% phi_x = zeros(num);
% xy = repmat(-floor(num/2):floor(num/2), num, 1);
% mask = xy.^2 + xy'.^2 < floor(num/2)^2;
% phi_x(mask) = 1; 
% phi_x(~mask) = -1;
% clear xy mask;

phi_x = -ones(num);
phi_x(floor(num/3):floor(2*num/3),:) = 1;


% 
% X=ones(num,1)*[-150:149];
% Y=[-150:149]'*ones(1,num);
% Z=X.^2+Y.^2;
% image=zeros([num num]);
% image(find(Z<=1000))=1;
% figure(2);imagesc(image);axis image;colormap(gray)
% IMAGE=1-image;

 %In order to use the reduced memory optimisation routine we need to switch
 %to a constrained optimisation.  So these are some 'dummy' constraints
 lb = -ones(size(phi_x));
 ub = ones(size(phi_x));
 
 options = optimset('fmincon');

% options = optimset(options,'MaxFunEvals',2*10^7,'MaxIter',2*10^7,'display','iter');

options = optimset(options,'MaxFunEvals',10^6,'MaxIter',2,'display','iter');
%options = optimset(options,'Algorithm','interior-point');
%options = optimset(options,'Hessian','on');

for ii = 1:10
    [phi_x] = fmincon(@ft_2d_complete, phi_x, [],[],[],[],lb,ub,[],options);
    figure; imagesc(phi_x); axis image;
    title(['Iteration: ' num2str(ii)]);
end
% [Fnew1]=fmincon(@ft_2d_secondterm_new,IMAGE,[],[],[],[],lb,ub,[],options);
%%
num = 127;
phi_x = zeros(num);
xy = repmat(-floor(num/2):floor(num/2), num, 1);
inner_mask = xy.^2 + xy'.^2 > floor(num/4)^2;

for ww = 1:floor(num/4)
    rr = ww + floor(num/4);
    outer_mask = xy.^2 + xy'.^2 < rr^2;
    
    mask = inner_mask & outer_mask;
    
    
    phi_x(mask) = 1; 
    phi_x(~mask) = -1;
    %figure; imagesc(phi_x); axis image;
    
    [gx gy] = gradient(phi_x);
    
    es_d2(ww) = CalculateES(gx, gy, 2); %#ok
    es_d4(ww) = CalculateES(gx, gy, 4); %#ok
    es_d8(ww) = CalculateES(gx, gy, 8); %#ok
    
    eo_d2(ww) = CalculateEO(phi_x, gx, gy); %#ok
    eo_d4(ww) = CalculateEO(phi_x, gx, gy); %#ok
    eo_d8(ww) = CalculateEO(phi_x, gx, gy); %#ok
end
%%
figure; hold on;
plot(1:floor(num/4), es_d2, 'r');
plot(1:floor(num/4), es_d4, 'g');
plot(1:floor(num/4), es_d8, 'b');

figure; hold on;
plot(1:floor(num/4), eo_d2, 'r');
plot(1:floor(num/4), eo_d4, 'g');
plot(1:floor(num/4), eo_d8, 'b');
%%
num = 127;
phi_x = zeros(num);
xy = repmat(-floor(num/2):floor(num/2), num, 1);
inner_mask = xy.^2 + xy'.^2 > floor(num/4)^2;

for ww = 1:floor(num/4)
    rr = ww + floor(num/4);
    outer_mask = xy.^2 + xy'.^2 < rr^2;
    
    mask = inner_mask & outer_mask;
    
    
    phi_x(mask) = 1; 
    phi_x(~mask) = -1;
    %figure; imagesc(phi_x); axis image;
    
    [gx gy] = gradient(phi_x);
    
    es_d2(ww) = CalculateES(gx, gy, 2); %#ok
    es_d4(ww) = CalculateES(gx, gy, 4); %#ok
    es_d8(ww) = CalculateES(gx, gy, 8); %#ok
    
    eo_d2(ww) = CalculateEO(phi_x, gx, gy); %#ok
    eo_d4(ww) = CalculateEO(phi_x, gx, gy); %#ok
    eo_d8(ww) = CalculateEO(phi_x, gx, gy); %#ok
end
%%
figure; hold on;
plot(1:floor(num/4), es_d2, 'r');
plot(1:floor(num/4), es_d4, 'g');
plot(1:floor(num/4), es_d8, 'b');

figure; hold on;
plot(1:floor(num/4), eo_d2, 'r');
plot(1:floor(num/4), eo_d4, 'g');
plot(1:floor(num/4), eo_d8, 'b');