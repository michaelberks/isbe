clc;

k = 5;

figure(1); clf; hold off; colormap(gray(256));
for i = 0:1
    if i == 0
        % Real
        imgroot = 'U:\projects\nailfold\capture\';
        imgroot = fullfile(imgroot, '2012_10_22');
            imgroot = fullfile(imgroot, 'Left.Digit4.x300\seq1\preprocessed\');
%             imgroot = fullfile(imgroot, 'Left.Digit4.x300\seq2\corrected\');
%         imgroot = fullfile(imgroot, '2013_02_21');
%             imgroot = fullfile(imgroot, 'Left.Digit4.x300\10_36_55\corrected\');

        imgroot = fullfile(imgroot, 'registered_g1d\');
        imgroot = fullfile(imgroot, 'cropped1');
%         imgroot = fullfile(imgroot, 'cropped2');
%         imgroot = fullfile(imgroot, 'cropped3');
        
        imgStack = load_image_stack(imgroot);
    else
        % Synthetic
        imgroot = 'U:\projects\nailfold\synthesis\';
        d = dir(fullfile(imgroot, '*T*'));
        imgroot = fullfile(imgroot, d(end).name);
%         imgroot = fullfile(imgroot, 'halfsize');
        
        imgStack = load_image_stack(imgroot);
        imgStack = imgStack(135:end-20, 23:end-22, :);
    end
    
    [m,n,nf] = size(imgStack);

    meanImage = mean(imgStack,3);
    stdImage = std(imgStack, [], 3);

    subplot(2,k,i*k+1);
        image(meanImage);
        axis('image','ij');
    subplot(2,k,i*k+2);
        image(uint8(stdImage*32));
        axis('image','ij');

    brightness = zeros(1,nf);
    contrast = zeros(1,nf);
    for j = 1:size(imgStack,3)
        imgStack(:,:,j) = (imgStack(:,:,j) - meanImage).^2;
        imgj = imgStack(:,:,j);
        brightness(j) = mean(imgj(:));
        contrast(j) = std(imgj(:));
    end

    subplot(2,k,i*k+3);
        hist(stdImage(:), 40);

    disp([min(meanImage(:)) max(meanImage(:))]);

    % Pick n pairs of randomly selected pixels
    N = numel(meanImage);
    n = 2000;
    inds = randperm(N);
    diffs = meanImage(inds(1:n)) - meanImage(inds(n+1:2*n));

    subplot(2,k,i*k+4);
%         hist(abs(diffs), 40);
        hist(meanImage(:), 256);

    subplot(2,k,i*k+5); hold on;
        plot(brightness,'b-');
        plot(contrast,'r-');
end