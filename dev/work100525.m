dot_im = zeros(256);
dot_im(128,128) = 1;
dt = dtwavexfm2(dot_im, 6);
dti = dt_to_full_image(dt);

for level = 1:6

    figure; 
    for ori = 1:6; 
        subplot(2,3,ori); imagesc(real(dti(:,:,ori,level))); colormap(gray(256)); 
        caxis([min(real(dt{level}(:))) max(real(dt{level}(:)))]); axis image; 
    end
    figure; 
    for ori = 1:6; 
        subplot(2,3,ori); imagesc(imag(dti(:,:,ori,level))); colormap(gray(256));
        caxis([min(imag(dt{level}(:))) max(imag(dt{level}(:)))]); axis image;
    end
end
%%
dot_im = zeros(256);
dot_im(128,128) = 1;
dt = dtwavexfm2b(dot_im, 6);
dti = dt_to_full_image(dt);

c_min = min([min(real(dt{6}(:))) min(imag(dt{6}(:)))]);
c_max = max([max(real(dt{6}(:))) max(imag(dt{6}(:)))]);
%Real parts
figure;

for ori = 1:6
    subplot(2,6,ori);  imagesc(real(dti(:,:,ori,6))); colormap(gray(256)); caxis([c_min c_max]); axis image;
    subplot(2,6,ori+6);  imagesc(imag(dti(:,:,ori,6))); colormap(gray(256)); caxis([c_min c_max]); axis image;
    
    %write_im_from_colormap(real(dti(:,:,ori,6)), ['k:\isbe\misc\dt_response_real_' num2str(ori) '.bmp'], gray(256));
    %write_im_from_colormap(imag(dti(:,:,ori,6)), ['k:\isbe\misc\dt_response_imag_' num2str(ori) '.bmp'], gray(256));
end
%%
dot_im = zeros(256);
dot_im(128,128) = 1;
dt = dtwavexfm2(dot_im, 6, 'near_sym_b', 'qshift_b');

c_min = min([min(real(dt{6}(:))) min(imag(dt{6}(:)))]);
c_max = max([max(real(dt{6}(:))) max(imag(dt{6}(:)))]);
%Real parts
figure;
subplot(2,6,1); imagesc(imresize(imag(dt{6}(:,:,4)),4,'bilinear')); caxis([c_min c_max]); axis image;
subplot(2,6,2); imagesc(imresize(real(dt{6}(:,:,2)),4,'bilinear')); caxis([c_min c_max]); axis image;
subplot(2,6,3); imagesc(imresize(imag(dt{6}(:,:,6)),4,'bilinear')); caxis([c_min c_max]); axis image;
subplot(2,6,4); imagesc(imresize(real(dt{6}(:,:,1)),4,'bilinear')); caxis([c_min c_max]); axis image;
subplot(2,6,5); imagesc(imresize(imag(dt{6}(:,:,5)),4,'bilinear')); caxis([c_min c_max]); axis image;
subplot(2,6,6); imagesc(imresize(real(dt{6}(:,:,3)),4,'bilinear')); caxis([c_min c_max]); axis image;

subplot(2,6,7); imagesc(imresize(real(dt{6}(:,:,4)),4,'bilinear')); caxis([c_min c_max]); axis image;
subplot(2,6,8); imagesc(imresize(imag(dt{6}(:,:,2)),4,'bilinear')); caxis([c_min c_max]); axis image;
subplot(2,6,9); imagesc(imresize(real(dt{6}(:,:,6)),4,'bilinear')); caxis([c_min c_max]); axis image;
subplot(2,6,10); imagesc(imresize(imag(dt{6}(:,:,1)),4,'bilinear')); caxis([c_min c_max]); axis image;
subplot(2,6,11); imagesc(imresize(real(dt{6}(:,:,5)),4,'bilinear')); caxis([c_min c_max]); axis image;
subplot(2,6,12); imagesc(imresize(imag(dt{6}(:,:,3)),4,'bilinear')); caxis([c_min c_max]); axis image;

%%
dot_im = zeros(256);
dot_im(128,128) = 1;
dt = dtwavexfm2(dot_im, 6, 'near_sym_b', 'qshift_d');

for level = 5:6

    figure; 
    for ori = 1:6; 
        subplot(2,3,ori); mesh(imresize(real(dt{level}(:,:,ori)),4,'bilinear')); 
        %zlims([min(real(dt{level}(:))) max(real(dt{level}(:)))]); 
    end
    figure; 
    for ori = 1:6; 
        subplot(2,3,ori); mesh(imresize(imag(dt{level}(:,:,ori)),4,'bilinear'));
        %zlims([min(imag(dt{level}(:))) max(imag(dt{level}(:)))]);
    end
end
%%
dot_im = zeros(256);
dot_im(128,128) = 1;
dt = dtwavexfm2(dot_im, 6, 'near_sym_b', 'qshift_d');

for level = 1:6

    figure; 
    for ori = 1:6; 
        subplot(2,3,ori); imagesc(imresize(abs(dt{level}(:,:,ori)),1,'bilinear')); %colormap(gray(256)); 
        caxis([0 max(abs(dt{level}(:)))]); axis image; 
    end
end