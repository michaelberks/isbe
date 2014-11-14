tag1 = imread('C:\isbe\density\IDtags\3A_10.bmp');
tag2 = imread('C:\isbe\density\IDtags\3A_6.bmp');
tag3 = imread('C:\isbe\density\IDtags\3A_2.bmp');
%%
figure; imagesc(tag1); axis image;
figure; imagesc(tag2); axis image;
figure; imagesc(tag3); axis image;

figure; hist(tag1(:), 0:255);
figure; hist(tag2(:), 0:255);
figure; hist(tag3(:), 0:255);
%%
tag1a = (double(tag1 - min(tag1(:)))) / double(max(tag1(:)) - min(tag1(:)));
tag2a = (double(tag2 - min(tag2(:)))) / double(max(tag2(:)) - min(tag2(:)));
tag3a = (double(tag3 - min(tag3(:)))) / double(max(tag3(:)) - min(tag3(:)));
figure; hist(tag1a(:), 100);
figure; hist(tag2a(:), 100);
figure; hist(tag3a(:), 100);
%%
figure; imagesc(tag1a); axis image; colormap(gray(256));
figure; imagesc(tag2a); axis image; colormap(gray(256));
figure; imagesc(tag3a); axis image; colormap(gray(256));
%%
tag1h = histeq(tag1,16);
tag2h = histeq(tag2,16);
tag3h = histeq(tag3,16);
figure; hist(tag1h(:), 0:255);
figure; hist(tag2h(:), 0:255);
figure; hist(tag3h(:), 0:255);
%%
figure; imagesc(histeq(tag1,2)); axis image; colormap(gray(256));
figure; imagesc(histeq(tag2,2)); axis image; colormap(gray(256));
figure; imagesc(histeq(tag3,2)); axis image; colormap(gray(256));
%%
imwrite(histeq(tag1,16), 'C:\isbe\density\IDtags\tag1h.bmp');
imwrite(histeq(tag2,16), 'C:\isbe\density\IDtags\tag2h.bmp');
imwrite(histeq(tag3,16), 'C:\isbe\density\IDtags\tag3h.bmp');

imwrite(tag1a, 'C:\isbe\density\IDtags\tag1a.bmp');
imwrite(tag2a, 'C:\isbe\density\IDtags\tag2a.bmp');
imwrite(tag3a, 'C:\isbe\density\IDtags\tag3a.bmp');

save('C:\isbe\density\IDtags\tag1a.mat', 'tag1a');
save('C:\isbe\density\IDtags\tag2a.mat', 'tag2a');
save('C:\isbe\density\IDtags\tag3a.mat', 'tag3a');
%%
for ii = 0:9
    digit = double(rgb2gray(imread(['C:\isbe\density\IDtags\' num2str(ii) '.bmp'])));
    digita = (digit - min(digit(:))) / (max(digit(:)) - min(digit(:)));
    digita(digita > 0.6) = 1;
%     figure;
%     subplot(1,2,1); imagesc(digit); axis image; colormap(gray(256));
%     subplot(1,2,2); imagesc(digita); axis image; colormap(gray(256));
    save(['C:\isbe\density\IDtags\' num2str(ii) '.mat'], 'digita');
end
%%
letters = 'ABCEGILMNOPRSTUY1234';
for ii = 1:length(letters)
    letter = double(rgb2gray(imread(['C:\isbe\density\IDtags\' letters(ii) '_big.bmp'])));
    letter = (letter - min(letter(:))) / (max(letter(:)) - min(letter(:)));
    figure;
    subplot(1,2,1); imagesc(letter); axis image; colormap(gray(256));
    letter(letter > 0.6) = 1;
    subplot(1,2,2); imagesc(letter); axis image; colormap(gray(256));
    save(['C:\isbe\density\IDtags\' letters(ii) '_big.mat'], 'letter');
end
%%
for tt = 1:3
    tag = u_load(['C:\isbe\density\IDtags\tag' num2str(tt) 'a.mat']);
    [rt ct] = size(tag);
    cc = zeros([rt ct 10]);
    for ii = 0:9
        dd = u_load(['C:\isbe\density\IDtags\' num2str(ii) '.mat']);
        [rd cd] = size(dd);
        rd2 = floor(rd/2);
        cd2 = floor(cd/2);
        c_temp = normxcorr2(dd, tag);
        cc(:,:,ii+1) = c_temp(rd2+(1:rt), cd2+(1:ct));
    %     figure;
    %     subplot(2,1,1); imagesc(tag1a); axis image;
    %     subplot(2,1,2); imagesc(cc); axis image;
    end
    % figure;
    % subplot(2,1,1); imagesc(tag1a); axis image;
    % subplot(2,1,2); imagesc(max(cc,[],3)); axis image;
    %
    [c_max c_id] = max(cc, [], 3);
    [maxima_pos, maxima_vals] = local_image_maxima(c_max, 10, [], 0.7);
    %
    figure; imagesc(tag); axis image; colormap(gray(256)); hold on;
    for ii = 1:size(maxima_pos,1)
        text(maxima_pos(ii,1), maxima_pos(ii,2), ...
            num2str(c_id(maxima_pos(ii,2), maxima_pos(ii,1))-1),...
            'color', 'r', 'fontsize', 16);
    end
end
%%
letters = 'ABCEGILMNOPRSTUY1234';
for tt = 1:3
    tag = u_load(['C:\isbe\density\IDtags\tag' num2str(tt) 'a.mat']);
    [rt ct] = size(tag);
    cc = zeros([rt ct 5]);
    
    for ii = 1:length(letters)
        dd = u_load(['C:\isbe\density\IDtags\' letters(ii) '_big.mat']);
        [rd cd] = size(dd);
        rd2 = floor(rd/2);
        cd2 = floor(cd/2);
        c_temp = normxcorr2(dd, tag);
        cc(:,:,ii) = c_temp(rd2+(1:rt), cd2+(1:ct));
    %     figure;
    %     subplot(2,1,1); imagesc(tag1a); axis image;
    %     subplot(2,1,2); imagesc(cc); axis image;
    end
%     figure;
%     subplot(2,1,1); imagesc(tag1a); axis image;
%     subplot(2,1,2); imagesc(max(cc,[],3)); axis image;
    %
    letter_mask = false(rt, ct);
    letter_mask(1:100,:) = 1;
    letter_mask(:,1100:end) = 1;
    [c_max c_id] = max(cc, [], 3);
    [maxima_pos] = local_image_maxima(c_max, 10, letter_mask, 0.7);
    %
    %letter_mask = true(rt, ct);
    figure; imagesc(tag); axis image; colormap(gray(256)); hold on;
    for ii = 1:size(maxima_pos,1)
        x = maxima_pos(ii,1);
        y = maxima_pos(ii,2);
        text(x-10, y, letters(c_id(y, x)),...
            'color', 'r', 'fontsize', 28);
        %letter_mask(y+(-40:40), x+(-20:20)) = 0;
    end
    
    for ii = 0:9
        dd = u_load(['C:\isbe\density\IDtags\' num2str(ii) '.mat']);
        [rd cd] = size(dd);
        rd2 = floor(rd/2);
        cd2 = floor(cd/2);
        c_temp = normxcorr2(dd, tag);
        cc(:,:,ii+1) = c_temp(rd2+(1:rt), cd2+(1:ct));
    end

    [c_max c_id] = max(cc, [], 3);
    [maxima_pos] = local_image_maxima(c_max, 10, ~letter_mask, 0.7);

    for ii = 1:size(maxima_pos,1)
        x = maxima_pos(ii,1);
        y = maxima_pos(ii,2);
        text(x-5, y, num2str(c_id(y, x)-1),...
            'color', 'r', 'fontsize', 16);
    end
end
%%
%Check function works
for tt = 1:3
    tag = u_load(['C:\isbe\density\IDtags\tag' num2str(tt) 'a.mat']);
    read_patient_id(tag, 'C:\isbe\density\IDtags\', 1);
end