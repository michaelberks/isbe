function [mask, pectoral, nipple] = convertPolarNipplePectLine(objLM, mask, lineprm, nipple, im_crop, fig)
    [N,M] = size(mask);
    if numel(lineprm) > 0
        [x1,y1,x2,y2] = polar2borderpoints([N M],lineprm(1),lineprm(2));
    end

    %mask(find(mask==0))=3; %#ok<FNDSB> %Change the label of the
    %background; SB - Commented by JM, no hard code!

    % Zero pad cropped mask
    mask = padarray(mask,[im_crop(1) im_crop(3)],objLM.labelBackground,'pre');
    mask = padarray(mask,[im_crop(2) im_crop(4)],objLM.labelBackground,'post');

    % Adjust points
    if numel(lineprm) > 0
        x1 = x1 + im_crop(1);
        x2 = x2 + im_crop(1);
        nipple(1) = nipple(1) + im_crop(1);

        y1 = y1 + im_crop(3);
        y2 = y2 + im_crop(3);
        nipple(2) = nipple(2) + im_crop(3);

        pectoral = [x1,y1;x2,y2];
    else
        pectoral = [];
    end

    % Figure
    if fig > 0
        figure;
%         subplot(1,3,1); imagesc(image); colormap('gray'); axis image; title('image');
%         if numel(lineprm) > 0
%             hold on; plot(nipple(2),nipple(1),'rs'); plot([y1,y2],[x1,x2]);
%             plot([y1,y2],[x1,x2],'rs'); hold off;
%         end
%         subplot(1,3,2); 
        imagesc(mask); colormap('gray'); axis image; title('mask');
        if numel(lineprm) > 0
            hold on; plot(nipple(2),nipple(1),'rs'); hold off;
            hold on; plot(nipple(2),nipple(1),'rs'); plot([y1,y2],[x1,x2]);
            plot([y1,y2],[x1,x2],'rs'); hold off;
        end
    end

end