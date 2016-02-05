function im_direction =  maskDirection(mask,labelBackground)
    tmpmask = mask;
    tmpmask(mask == labelBackground) = 0;
    
    numberColumns = size(tmpmask,2);
    iml = tmpmask(:,1:round(numberColumns/2));
    imr = tmpmask(:,round(numberColumns/2):end);
    
    if ( mean(imr(:)) >  mean(iml(:)) ) 
        im_direction = 'right';
    else
        im_direction = 'left';
    end

end