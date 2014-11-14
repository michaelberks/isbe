bg = double(imread(['C:\isbe\dev\background\images\normal512\o04_004RCC_1024_5357_4579.bmp']));
%%
[g_pyr] = mb_g_pyramid(bg, 5);
%%
for ii = 1:5
    figure; imagesc(g_pyr{ii}); axis image; colormap(gray(256));
end
%%
[linop_pyr ori_pyr] = line_operator_pyr(bg, 12, 5);
%%
for ii = 1:5
    figure; imagesc(linop_pyr{ii}); axis image; colormap(gray(256));
end
%%
[linop_pyr_full] = line_operator_pyr_full(bg, 12, 5);
for ii = 1:5
    figure; 
    for ori = 1:12
        subplot(3,4,ori);
        imagesc(linop_pyr_full{ii}(:,:,ori)); axis image; colormap(gray(256));
    end
end
%%
linop_response = line_operator_octave(bg, 12, 5, 'degrees', 1);
figure; imagesc(linop_response); axis image; colormap(gray(256));
%%
for ii = 1:100
    test_image = u_load(['M:\chen\data\testimage_contrast0to8_multibars_sin\image' num2str(ii) '.mat']);
    load(['M:\chen\data\testimage_contrast0to8_multibars_sin\label' num2str(ii) '.mat']);
    label_centre = multibar_labelcentre;
    
    save(['C:\isbe\dev\classification\data\testimage_contrast0to8_multibars_sin\test_image' zerostr(ii,3) '.mat'], 'test_image', 'label_centre');
end
%%
for ii = 1:100
    test_image = u_load(['M:\chen\data\testimage_contrast0to8_multibars_sin\image' num2str(ii) '.mat']);
    load(['M:\chen\data\testimage_contrast0to8_multibars_sin\label' num2str(ii) '.mat']);
    label_centre = multibar_labelcentre;
    label = multibar_label;
    
    save(['C:\isbe\dev\classification\data\testimage_contrast0to8_multibars_sin\test_image' zerostr(ii,3) '.mat'], 'test_image', 'label_centre', 'label');
end
%%
for ii = 1:100
    probability_image = ...
        u_load(['C:\isbe\dev\classification\data\testimage_contrast0to8_multibars_sin\probability_images\182269\probability_image' num2str(ii) '.mat']);
    probability_image = 1 - probability_image;
    save(['C:\isbe\dev\classification\data\testimage_contrast0to8_multibars_sin\probability_images\182269\probability_image' zerostr(ii,3) '.mat'],...
        'probability_image');
    delete(['C:\isbe\dev\classification\data\testimage_contrast0to8_multibars_sin\probability_images\182269\probability_image' num2str(ii) '.mat']);
        
    probability_image = ...
        u_load(['C:\isbe\dev\classification\data\testimage_contrast0to8_multibars_sin\probability_images\182270\probability_image' num2str(ii) '.mat']);
    probability_image = 1 - probability_image;
    save(['C:\isbe\dev\classification\data\testimage_contrast0to8_multibars_sin\probability_images\182270\probability_image' zerostr(ii,3) '.mat'],...
        'probability_image');
    delete(['C:\isbe\dev\classification\data\testimage_contrast0to8_multibars_sin\probability_images\182270\probability_image' num2str(ii) '.mat']);
    
    probability_image = ...
        u_load(['C:\isbe\dev\classification\data\testimage_contrast0to8_multibars_sin\probability_images\182272\probability_image' num2str(ii) '.mat']);
    probability_image = 1 - probability_image;
    save(['C:\isbe\dev\classification\data\testimage_contrast0to8_multibars_sin\probability_images\182272\probability_image' zerostr(ii,3) '.mat'],...
        'probability_image');
    delete(['C:\isbe\dev\classification\data\testimage_contrast0to8_multibars_sin\probability_images\182272\probability_image' num2str(ii) '.mat']);

end
%%
for ii = 1:100
    probability_image = ...
        u_load(['C:\isbe\dev\classification\data\testimage_contrast0to8_multibars_sin\probability_images\182267\probability_image' zerostr(ii,3) '.mat']);
    probability_image = 1 - probability_image;
    save(['C:\isbe\dev\classification\data\testimage_contrast0to8_multibars_sin\probability_images\182267\probability_image' zerostr(ii,3) '.mat'],...
        'probability_image');
end