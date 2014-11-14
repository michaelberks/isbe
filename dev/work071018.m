idx = randsample(101, 20);

for ii = 1:20
    jj = idx(ii);
    mass_image = zeros(mass_model.mean_row, mass_model.mean_col);
    mass_image(sub2ind(size(mass_image), mass_model.mean_shape_pl(:,1),...
        mass_model.mean_shape_pl(:,2))) = mass_model.X_tex(jj,:);
    figure; imagesc(mass_image); colormap(gray(256)); axis image;
end

%%
for ii = 1:10
    mass_image = new_masses(ii).mass_ROI;
    padding = floor((size(bg) - size(mass_image))/2);
    if all(padding > 0)
        nm = bg;
        nm(padding(1):padding(1)+size(mass_image,1)-1, padding(2):padding(2)+size(mass_image,2)-1) = nm(padding(1):padding(1)+size(mass_image,1)-1, padding(2):padding(2)+size(mass_image,2)-1) + mass_image;
        figure; imagesc(nm); colormap(gray(256)); axis image;
    end
end
    