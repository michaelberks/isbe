function view_subject(subject_id, varargin)

args = u_packargs(varargin,... % the user's input
    0, ... % non-strict mode
    'image_dir', 'C:\isbe\nailfold\data\rsa_study\images\',...
    'im_fmt', 'png',...
    'visit', 1);

load C:\isbe\nailfold\data\rsa_study\image_id_data.mat

%i_im = find(strcmp(image_id_data.im_names, '10598c'));
%pp = image_id_data.people_id(i_im);
%vi = image_id_data.visit(i_im);

for pp = subject_id(:)'


    idx = find((image_id_data.people_id == pp) & (image_id_data.visit == args.visit) );% & (image_id_data.digit >= 2)

    
    figure(...
        'name', ['Subject ' num2str(pp), ' from group ' image_id_data.category{idx(1)} ', visit ' num2str(args.visit)]);
    
    for i_im = 1:length(idx)
        im_name = image_id_data.im_names{idx(i_im)};
        hand = image_id_data.hand{idx(i_im)};
        digit = image_id_data.digit(idx(i_im));
        
        if exist([args.image_dir im_name '.' args.im_fmt], 'file')
            im = imread(['C:\isbe\nailfold\data\rsa_study\images\' im_name '.png']);
            im = im(:,:,1);
            mask = make_nailfold_mosaic_mask(im);
        else
            continue
        end
        
%         if exist(['C:\isbe\nailfold\data\rsa_study\test_half\images\' im_name '.mat'], 'file')
%             im = u_load(['C:\isbe\nailfold\data\rsa_study\test_half\images\' im_name '.mat']);
%             mask = u_load(['C:\isbe\nailfold\data\rsa_study\test_half\fov_masks\' im_name '_f_mask.mat']);
%         elseif exist(['C:\isbe\nailfold\data\rsa_study\final_test\images\' im_name '.mat'], 'file')
%             im = u_load(['C:\isbe\nailfold\data\rsa_study\final_test\images\' im_name '.mat']);
%             mask = u_load(['C:\isbe\nailfold\data\rsa_study\final_test\fov_masks\' im_name '_f_mask.mat']);
%         elseif exist(['C:\isbe\nailfold\data\rsa_study\images\' im_name '.png'], 'file')
%             im = imread(['C:\isbe\nailfold\data\rsa_study\images\' im_name '.png']);
%             im = im(:,:,1);
%             mask = true(size(im));
%         elseif exist(['Q:\nailfold\data\rsa_study\images\caps\' im_name '.bmp'], 'file')
%             im = imread(['Q:\nailfold\data\rsa_study\images\caps\' im_name '.bmp']);
%             im = im(:,:,1);
%             mask = true(size(im));
%         else
%             continue;
%         end

        plot_r = 5 - digit;
        if hand == 'L'
            plot_c = 0;
        else
           plot_c = 0.5;
        end
        
        axes('position', [plot_c 0.2*plot_r 0.5 0.2]);
        
        imgray(im);
        %title([im_name ' H: ' hand ', D: ' num2str(digit) ', V: ' num2str(args.visit)]);
        text(10, 20, im_name, 'color', 'r');
        g_min = min(im(mask)) + 5;
        g_max = max(im(mask)) - 5;
        caxis([g_min g_max]);
        axis off;
    end
end
