%%
f1 = figure('WindowStyle', 'docked'); image(dense); colormap(gray(256)); axis('image'); hold('on');
f2 = figure('WindowStyle', 'docked'); image(mixed1); colormap(gray(256)); axis('image'); hold('on');
f3 = figure('WindowStyle', 'docked'); image(mixed2); colormap(gray(256)); axis('image'); hold('on');
f4 = figure('WindowStyle', 'docked'); image(fatty); colormap(gray(256)); axis('image'); hold('on');
%%
mass16(16).outline_500 = [];
for ii = 1:16
   mass_ROI = mass16(ii).mass_ROI;
   mass_outline = mass16(ii).mass_outline;
   if ii == 8; mass_outline = mass_outline * 0.9; end
   c_x = min(mass_outline(:,1)) +...
       (max(mass_outline(:,1)) - min(mass_outline(:,1)))/2;
   c_y = min(mass_outline(:,2)) +...
       (max(mass_outline(:,2)) - min(mass_outline(:,2)))/2;
   
   mass_outline(:,1) = mass_outline(:,1) + 500 - c_x;
   mass_outline(:,2) = mass_outline(:,2) + 500 - c_y;
   
   figure(f1);
   plot(mass_outline(:,1), mass_outline(:,2));
   figure(f2);
   plot(mass_outline(:,1), mass_outline(:,2));
   figure(f3);
   plot(mass_outline(:,1), mass_outline(:,2));
   figure(f4);
   plot(mass_outline(:,1), mass_outline(:,2));
   mass16(ii).outline_500 = mass_outline;
end
%%

inner_area(16) = 0;
outer_area(16) = 0;
bg_var(16) = 0;
mass_var(16) = 0;
bg_mean(16) = 0;
mass_mean(16) = 0;

for ii = 1:16
   mass_ROI = mass16(ii).mass_ROI;
   mass_outline = mass16(ii).mass_outline;
   
   mass_bw = roipoly(mass_ROI, mass_outline(:,1), mass_outline(:,2));
   inner_area(ii) = sum(mass_bw(:));
   for jj = 1:20; mass_bw = imdilate(mass_bw, strel('disk', 1)); end
   
   mass_pl = regionprops(bwlabel(mass_bw, 4), 'PixelIdxList');
   p_list = sort(mass_pl.PixelIdxList); clear mass_pl;
   
   mass_mean(ii) = mean(mass_ROI(p_list));
   bg_mean(ii) = mean(mass_ROI(setdiff(1:end, p_list)));
   
   t1 = imfilter(double(mass_ROI), ones(9)/81, 'symmetric');
   t2 = imfilter(double(mass_ROI).^2, ones(9)/81, 'symmetric');
   local_var = t2 - (t1.*t1);
   
   mass_var(ii) = mean(local_var(p_list));
   bg_var(ii) = mean(local_var(setdiff(1:end, p_list)));
end

%%
inner_area(16) = 0;
for ii = 1:16
   mass_outline = mass_to_sub(ii).mass_outline;
   
   mass_bw = roipoly(1000,1000, mass_outline(:,1), mass_outline(:,2));
   inner_area(ii) = sum(mass_bw(:));
end
clear mass_bw mass_outline
   