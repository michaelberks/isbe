for ii = 1:16
    
    frames(:,:,ii) = imread(['E:\nailfold\visit1\Lhand\Digit4\X300\Scene1\2V1LD4X3S1F' num2str(ii) '.bmp']);
end
%%
mask = f1 < 100;
figure; imagesc(mask); axis image;
f1(mask) = mean(f1(~mask));
figure; imagesc(f1); axis image; colormap(gray(256)); colorbar
%%

for win_size = [11 33 55]
    snr = mb_snr(f1, win_size);
    display(mean(snr(:)));
end