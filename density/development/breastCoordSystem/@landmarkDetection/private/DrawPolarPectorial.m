function DrawPolarPectorial(im, mask, lineprm, nipple)
    figure;
    subplot(1,2,1); imshow(im,[]); title('image');
    if numel(lineprm) > 0
        DrawLines_Polar(size(im), lineprm);
        hold on; plot(nipple(2),nipple(1),'rs'); hold off;
    end
    subplot(1,2,2); imshow(mask,[]); title('mask');
    if numel(lineprm) > 0
        DrawLines_Polar(size(im), lineprm);
        hold on; plot(nipple(2),nipple(1),'rs'); hold off;
    end
end