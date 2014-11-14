function display_modes_spicule(spicule_AM)    
%%%%%%%%%%%%%%%%%
%
% Function probably obsolete but may have some useful salvagable code!
%
%%%%%%%%%%%%%%%%%
    mean_s  = spicule_AM.mean_s;
    P_s     = spicule_AM.P_s;
    B_s     = spicule_AM.B_s;
    L_s     = spicule_AM.L_s;

    mean_c  = spicule_AM.mean_c;
    P_c     = spicule_AM.P_c;
    B_c     = spicule_AM.B_c;
    L_c     = spicule_AM.L_c;
            
    shape1_p = mean_s + (2*P_s(:,1)*sqrt(L_s(1)))';
    shape1_m = mean_s + (-2*P_s(:,1)*sqrt(L_s(1)))';

    shape2_p = mean_s + (2*P_s(:,2)*sqrt(L_s(2)))';
    shape2_m = mean_s + (-2*P_s(:,2)*sqrt(L_s(2)))';
    
    shape3_p = mean_s + (2*P_s(:,3)*sqrt(L_s(3)))';
    shape3_m = mean_s + (-2*P_s(:,3)*sqrt(L_s(3)))';
    
    c1_p = ((mean_c + (2*P_c(:,1)*sqrt(L_c(1)))') * P_c)';
    [c_shape1_p c_ROI1_p] = make_spicule(c1_p);
    
    c1_m = ((mean_c + (-2*P_c(:,1)*sqrt(L_c(1)))') * P_c)';
    [c_shape1_m c_ROI1_m] = make_spicule(c1_m);
    
    c2_p = ((mean_c + (2*P_c(:,2)*sqrt(L_c(2)))') * P_c)';
    [c_shape2_p c_ROI2_p] = make_spicule(c2_p);
    
    c2_m = ((mean_c + (-2*P_c(:,2)*sqrt(L_c(2)))') * P_c)';
    [c_shape2_m c_ROI2_m] = make_spicule(c2_m);
    
    c4_p = ((mean_c + (2*P_c(:,4)*sqrt(L_c(4)))') * P_c)';
    [c_shape4_p c_ROI4_p] = make_spicule(c4_p);
    
    c4_m = ((mean_c + (-2*P_c(:,4)*sqrt(L_c(4)))') * P_c)';
    [c_shape4_m c_ROI4_m] = make_spicule(c4_m);
    
    c3_p = ((mean_c + (2*P_c(:,3)*sqrt(L_c(3)))') * P_c)';
    [c_shape3_p c_ROI3_p] = make_spicule(c3_p);
    
    c3_m = ((mean_c + (-2*P_c(:,3)*sqrt(L_c(3)))') * P_c)';
    [c_shape3_m c_ROI3_m] = make_spicule(c3_m);
    
    [c_mean c_ROI] = make_spicule(zeros(50,1));
    
    figure,
    subplot(3,3,1)
    plot(shape1_m(1:end/2), shape1_m(1+end/2:end)); axis equal; 
    subplot(3,3,2)
    plot(mean_s(1:end/2), mean_s(1+end/2:end)); axis equal;
    subplot(3,3,3)
    plot(shape1_p(1:end/2), shape1_p(1+end/2:end)); axis equal; 

    subplot(3,3,4)
    plot(shape2_m(1:end/2), shape2_m(1+end/2:end)); axis equal; 
    subplot(3,3,5)
    plot(mean_s(1:end/2), mean_s(1+end/2:end)); axis equal;
    subplot(3,3,6)
    plot(shape2_p(1:end/2), shape2_p(1+end/2:end)); axis equal;
    
    subplot(3,3,7)
    plot(shape3_m(1:end/2), shape3_m(1+end/2:end)); axis equal; 
    subplot(3,3,8)
    plot(mean_s(1:end/2), mean_s(1+end/2:end)); axis equal;
    subplot(3,3,9)
    plot(shape3_p(1:end/2), shape3_p(1+end/2:end)); axis equal;

    figure,
    subplot(3,3,1)
    imagesc(c_ROI1_m); colormap gray; hold on; axis image;
    subplot(3,3,2)
    imagesc(c_ROI); colormap gray; hold on; axis image;
    subplot(3,3,3)
    imagesc(c_ROI1_p); colormap gray; hold on; axis image;

    subplot(3,3,4)
    imagesc(c_ROI2_m); colormap gray; hold on; axis image; 
    subplot(3,3,5)
    imagesc(c_ROI); colormap gray; hold on; axis image;
    subplot(3,3,6)
    imagesc(c_ROI2_p); colormap gray; hold on; axis image;
    
    subplot(3,3,7)
    imagesc(c_ROI3_m); colormap gray; hold on; axis image; 
    subplot(3,3,8)
    imagesc(c_ROI); colormap gray; hold on; axis image;
    subplot(3,3,9)
    imagesc(c_ROI3_p); colormap gray; hold on; axis image;
    
    subplot(4,3,10)
    imagesc(c_ROI4_m); colormap gray; hold on; axis image; 
    subplot(4,3,11)
    imagesc(c_ROI); colormap gray; hold on; axis image;
    subplot(4,3,12)
    imagesc(c_ROI4_p); colormap gray; hold on; axis image;
    
    figure('Position', [100,100, 618, 758]);
    aviobj = avifile('K:\isbe\project\misc\ss07\spic_modes.avi','fps', 10);
	for ii = 1:20
        k = -2 + 4*ii/20;
        
        subplot(3,2,1);
        shape_ii = mean_s + (k*P_s(:,1)*sqrt(L_s(1)))';
        plot(shape_ii(1:end/2), shape_ii(1+end/2:end)); axis equal;
        title('1st Mode of Shape');
        
        subplot(3,2,2);
        c_ii = ((mean_c + (k*P_c(:,1)*sqrt(L_c(1)))') * P_c)';
        [c_shape_ii c_ROI_ii] = make_spicule(c_ii);
        imagesc(c_ROI_ii); colormap gray; hold on; axis image;
        title('1st Combined Appearance Mode');
        
        subplot(3,2,3);
        shape_ii = mean_s + (k*P_s(:,2)*sqrt(L_s(2)))';
        plot(shape_ii(1:end/2), shape_ii(1+end/2:end)); axis equal;
        title('2nd Mode of Shape');
        
        subplot(3,2,4);
        c_ii = ((mean_c + (k*P_c(:,2)*sqrt(L_c(2)))') * P_c)';
        [c_shape_ii c_ROI_ii] = make_spicule(c_ii);
        imagesc(c_ROI_ii); colormap gray; hold on; axis image;
        title('2nd Combined Appearance Mode');
        
        subplot(3,2,5);
        shape_ii = mean_s + (k*P_s(:,3)*sqrt(L_s(3)))';
        plot(shape_ii(1:end/2), shape_ii(1+end/2:end)); axis equal;
        title('3rd Mode of Shape');
        
        subplot(3,2,6);
        c_ii = ((mean_c + (k*P_c(:,3)*sqrt(L_c(3)))') * P_c)';
        [c_shape_ii c_ROI_ii] = make_spicule(c_ii);
        imagesc(c_ROI_ii); colormap gray; hold on; axis image;
        title('3rd Combined Appearance Mode');
        
        frame = getframe(gcf);
        aviobj = addframe(aviobj,frame);
    end
    for ii = 1:20
        k = 2 - 4*ii/20;
        
        subplot(3,2,1);
        shape_ii = mean_s + (k*P_s(:,1)*sqrt(L_s(1)))';
        plot(shape_ii(1:end/2), shape_ii(1+end/2:end)); axis equal;
        title('1st Mode of Shape');
        
        subplot(3,2,2);
        c_ii = ((mean_c + (k*P_c(:,1)*sqrt(L_c(1)))') * P_c)';
        [c_shape_ii c_ROI_ii] = make_spicule(c_ii);
        imagesc(c_ROI_ii); colormap gray; hold on; axis image;
        title('1st Combined Appearance Mode');
        
        subplot(3,2,3);
        shape_ii = mean_s + (k*P_s(:,2)*sqrt(L_s(2)))';
        plot(shape_ii(1:end/2), shape_ii(1+end/2:end)); axis equal;
        title('2nd Mode of Shape');
        
        subplot(3,2,4);
        c_ii = ((mean_c + (k*P_c(:,2)*sqrt(L_c(2)))') * P_c)';
        [c_shape_ii c_ROI_ii] = make_spicule(c_ii);
        imagesc(c_ROI_ii); colormap gray; hold on; axis image;
        title('2nd Combined Appearance Mode');
        
        subplot(3,2,5);
        shape_ii = mean_s + (k*P_s(:,3)*sqrt(L_s(3)))';
        plot(shape_ii(1:end/2), shape_ii(1+end/2:end)); axis equal;
        title('3rd Mode of Shape');
        
        subplot(3,2,6);
        c_ii = ((mean_c + (k*P_c(:,3)*sqrt(L_c(3)))') * P_c)';
        [c_shape_ii c_ROI_ii] = make_spicule(c_ii);
        imagesc(c_ROI_ii); colormap gray; hold on; axis image;
        title('3rd Combined Appearance Mode');
        
        frame = getframe(gcf);
        aviobj = addframe(aviobj,frame);
    end
    aviobj = close(aviobj); clear aviobj;
     
    function [new_shape, new_ROI] = make_spicule(new_c)
        
        mean_w  = spicule_AM.mean_w;
        P_w     = spicule_AM.P_w;
        B_w     = spicule_AM.B_w;
        L_w     = spicule_AM.L_w;
        mean_b  = spicule_AM.mean_b;
        P_b     = spicule_AM.P_b;
        B_b     = spicule_AM.B_b;
        L_b     = spicule_AM.L_b;
        mean_p  = spicule_AM.mean_p;
        P_p     = spicule_AM.P_p;
        B_p     = spicule_AM.B_p;
        L_p     = spicule_AM.L_p;
        mean_l  = spicule_AM.mean_l;
        P_l     = spicule_AM.P_l;
        B_l     = spicule_AM.B_l;
        L_l     = spicule_AM.L_l;

        W_s     = spicule_AM.W_s;
        W_w     = spicule_AM.W_w;
        W_b     = spicule_AM.W_b;
        W_p     = spicule_AM.W_p;
        W_l     = spicule_AM.W_l;

        m = length(mean_s) / 2;

        k_s = length(L_s);
        k_w = length(L_w);
        k_b = length(L_b);
        k_p = length(L_p);

        Q_s = P_c(1:k_s,:); 
        Q_w = P_c(k_s+1:k_s + k_w,:);
        Q_b = P_c(k_s+k_w+1:k_s+k_w+k_b,:);
        Q_p = P_c(k_s+k_w+k_b+1:end-1,:);
        Q_l = P_c(end, :);

        new_s = mean_s + (P_s*Q_s*new_c)' / W_s;
        new_w = mean_w + (P_w*Q_w*new_c)' / W_w;
        new_b = mean_b + (P_b*Q_b*new_c)' / W_b;
        new_p = mean_p + (P_p*Q_p*new_c)' / W_p;
        new_l = mean_l + (P_l*Q_l*new_c)' / W_l;

        new_b(new_b < 1) = 1;
        new_w(new_w < 1) = 1;
        
        new_shape(:,1) = new_s(1:m);
        new_shape(:,2) = new_s(m+1:2*m);
        new_shape = new_shape*new_l / sqrt(sum((new_shape(end,:) - new_shape(1,:)).^2));
        
        %calculate length of new spicule
        V = diff(new_shape);
        N = [0; cumsum(sqrt(V(:,1).^2 + V(:,2).^2))];

        %create landmark points spaced along pixel length for new shape and
        %profile
        n_pts = round(N(end) / 5);%
        landmarks = interp1(N, new_shape, linspace(0, N(end), n_pts));

        % Reconstruct profile appearance vectors
        profile_x = [1:N(end)+1];
        prows = 81; pcols = length(profile_x);

        landmarks_xp = linspace(1, floor(N(end) + 1), n_pts);
        landmarks_yp = repmat(41, 1, n_pts);
        landmarks_w = ...
            interp1(linspace(1, pcols, length(new_w)), new_w, landmarks_xp);    
        upper_yp = landmarks_yp + landmarks_w;
        lower_yp = landmarks_yp - landmarks_w;

        %calculate widths at landmarks along profile from new_w
        w_vector = interp1(linspace(1, pcols, length(new_w)), new_w, 1:pcols);
        w_matrix = repmat(w_vector, prows, 1);

        %calculate brightness along profile from new_b
        b_vector = interp1(linspace(1, pcols, length(new_b)), new_b, 1:pcols);
        b_matrix = repmat(b_vector, prows, 1);

        offsets = repmat([-40:40]', 1, pcols);
        offsets(abs(offsets) >= w_matrix) = w_matrix(abs(offsets) >= w_matrix);

        new_profile = b_matrix .* sqrt(1 - (offsets.^2 ./ w_matrix.^2));

        %create upper and lower border landmarks for spicule shape
        [fx, fy] = gradient(landmarks);
        fy = fy ./ [sqrt(sum((fy.^2)')'), sqrt(sum((fy.^2)')')];

        landmarks_x = landmarks(:,1)'; 
        landmarks_y = landmarks(:,2)'; clear landmarks;


        upper_x = landmarks_x - landmarks_w.*fy(:,2)';
        upper_y = landmarks_y + landmarks_w.*fy(:,1)';
        lower_x = landmarks_x + landmarks_w.*fy(:,2)';
        lower_y = landmarks_y - landmarks_w.*fy(:,1)';

        %translate xy-coords to fit in box
        xbuff = 10 - min([upper_x, lower_x]);
        upper_x = upper_x + xbuff;
        lower_x = lower_x + xbuff;
        landmarks_x = landmarks_x + xbuff;
        new_shape(:,1) = new_shape(:,1) + xbuff;

        ybuff = 10 - min([upper_y, lower_y]); %NB: if spicule "hooks" then min could be in lower or upper
        upper_y = upper_y + ybuff;
        lower_y = lower_y + ybuff;
        landmarks_y = landmarks_y + ybuff;
        new_shape(:,2) = new_shape(:,2) + ybuff;

        %create box big enough for new spicule to fit into
        new_r = round(max([upper_y, lower_y])) + 10;
        new_c = round(max([upper_x, lower_x])) + 10;
        new_ROI = zeros(new_r, new_c);

        %create pixel list for new spicule, mask outline is upper border
        %landmarks from 1 to end, then lower landmarks from end to 1
        new_bw = roipoly(new_r, new_c, [upper_x fliplr(lower_x)],...
            [upper_y fliplr(lower_y)]);

        rp = regionprops(bwlabel(new_bw, 4), 'PixelList');

        new_shape_pl = rp.PixelList; clear rp;

        % Compute TPS warp to map from mean to new shape
        %%%%
        %Define source points for TPS - as row vectors
        s_x = [upper_x landmarks_x lower_x];
        s_y = [upper_y landmarks_y lower_y];

        %Define points to be interpolated by TPS - as row vectors
        i_x = new_shape_pl(:,1)';
        i_y = new_shape_pl(:,2)';

        tps_L_inv = tps_weights(s_x, s_y);

        %Define displacement to target points
        z_x = [landmarks_xp landmarks_xp landmarks_xp];
        z_y = [upper_yp landmarks_yp lower_yp];

        %Compute displacement of interpolated points        
        f_x = tps_warp(tps_L_inv, s_x, s_y, z_x, i_x, i_y);
        f_y = tps_warp(tps_L_inv, s_x, s_y, z_y, i_x, i_y);

        new_tex = interp2(double(new_profile), f_x, f_y);            
        new_ROI(sub2ind([new_r, new_c], i_y, i_x)) = uint8(new_tex);
    end