% probability_dir = 'M:\chen\data\line_detection_mammo\results\DTCWT_rf_fulltrees_W3L5_200000\';
% angle_dir = 'M:\chen\data\line_detection_mammo\results\regtree_W3L5_chen\';
% save_path = 'M:\chen\data\line_detection_mammo\results\probability_orientation\';

probability_dir = 'M:\chen\data\testimage_contrast1to8_exprnd_sin\synthetic\probability\';
angle_dir = 'M:\chen\data\testimage_contrast1to8_exprnd_sin\synthetic\regression\';
save_path = 'M:\chen\data\testimage_contrast1to8_exprnd_sin\synthetic\probability_orientation\';

if ~exist(save_path, 'dir')
    mkdir(save_path);
end

flag=0; %[1 0];
% for fl=1:2
    distflag = flag(fl);
    for m =2:5
        roi_prob = u_load([probability_dir, 'probability_image', zerostr(m,3), '.mat']);
        roi_ori = u_load([angle_dir, 'probability_image', zerostr(m,3), '.mat']);
        [ROW, COL] = size(roi_prob);
        im_out=zeros(ROW, COL);
        for rw = 1:ROW
            for cl = 1:COL
                xy=[];
                K = -tan(roi_ori(rw, cl)*pi/180);
                xy(1, 1) = 1;
                xy(1, 2) = rw + K*(1-cl);
                xy(2, 1) = cl + (1-rw)/K;
                xy(2, 2) = 1;
                xy(3, 1) = COL;
                xy(3, 2) = rw + K*(COL -cl);
                xy(4, 1) = cl + (ROW - rw)/K;
                xy(4, 2) = ROW;
                sort_xy = sortrows(xy);
                mycoordinates = round(sort_xy(2:3, :));
                [myline,mycoords,outmat,rr,cc] = bresenham(roi_prob,mycoordinates, 0);
                index = round(length(myline)/2);
                IND=sub2ind([ROW, COL], rr, cc);
                IND0=sub2ind([ROW, COL], rw, cl);
                myline(:,:) = roi_prob(IND0);
                if distflag
                    index = round(length(myline)/2);
                    IND=sub2ind([ROW, COL], rr, cc);
                    IND0=sub2ind([ROW, COL], rw, cl);
                    index = find(IND==IND0);
                    dist_decay = 1./[meshgrid(1:length(myline), 1)-index].^2;
                    dist_decay(isinf(dist_decay))=1;
                    myline = myline.* dist_decay;
                end
                im_out(IND) = im_out(IND)+ myline;
                %         hold on; plot(cl, rw, 'r+');
            end
        end
        disp(m);
        if distflag
            save([save_path, 'image_dist', zerostr(m, 3), '.mat'], 'im_out');
        else
            save([save_path, 'image', zerostr(m, 3), '.mat'], 'im_out');
        end
    end
% end