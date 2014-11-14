line_width = 2;
line_contrast = 4;

test_im = zeros(64);
for i_angle = 0:15:165

    %Generate line
    [line, label, label_centre] =...
        create_ellipse_bar(line_width/2, line_contrast, i_angle, 64, 64, 32, 32);
    
    test_im = max(test_im, line);
end
    
% perform the decomposition
[pyr, pind] = mb_buildSFpyr(test_im, 3, 5);

% now convert the pyramids to the CJR format
pyr = mb_change_pyramid_form(pyr, pind);

for i_level = 1:5
    figure;
    if i_level == 1 || i_level == 5;
        imgray(pyr{i_level,1});
    else
        for i_ori = 1:6
            subplot(2,3,i_ori);
            imgray(pyr{i_level,i_ori});
        end
    end
end
        