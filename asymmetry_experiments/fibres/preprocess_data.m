%%
for data = {'training\', 'test\'};
    mkdir([asymmetryroot 'data\fibre\' data{1} 'images\']);
    mkdir([asymmetryroot 'data\fibre\' data{1} 'fibre_masks\']);
    mkdir([asymmetryroot 'data\fibre\' data{1} 'red_fibre_masks\']);
    mkdir([asymmetryroot 'data\fibre\' data{1} 'blue_fibre_masks\']);
    mkdir([asymmetryroot 'data\fibre\' data{1} 'junction_fibre_masks\']);
    mkdir([asymmetryroot 'data\fibre\' data{1} 'fov_masks\']);
    mkdir([asymmetryroot 'data\fibre\' data{1} 'orientation_maps\']);
    mkdir([asymmetryroot 'data\fibre\' data{1} 'fibre_masks_dilated\']);
end
%%
fibre_list = dir([asymmetryroot 'data\fibre\CCMdatabase_Analysed\*.bmp']);
num_ims = length(fibre_list);
%
win_size = 11;
for ii = 201:300

    %Load in the original fibre image and the annotated copy
    orig_im = imread([asymmetryroot 'data\fibre\CCMdatabase\' fibre_list(ii).name]);
    anno_im = imread([asymmetryroot 'data\fibre\CCMdatabase_Analysed\' fibre_list(ii).name]);
    
    %Flatten the original image
    if size(orig_im,3) == 3
        orig_im = rgb2gray(orig_im);
    end
    
    %Extact the red blue and green annotations
    r_map = (anno_im(:,:,1) > 0.9) & (anno_im(:,:,2) < 0.1) & (anno_im(:,:,3) < 0.1);
    g_map = (anno_im(:,:,1) < 0.1) & (anno_im(:,:,2) > 0.9) & (anno_im(:,:,3) < 0.1);
    b_map = (anno_im(:,:,1) < 0.1) & (anno_im(:,:,2) < 0.1) & (anno_im(:,:,3) > 0.9);
    
    %Convert into masks of the junctions and fibres
    junction_mask = bwmorph(g_map, 'thin', inf);
    fibre_mask = r_map | b_map | junction_mask;
    red_fibre_mask = r_map | junction_mask;
    blue_fibre_mask = b_map | junction_mask;
    fibre_mask_d3 = bwdist(fibre_mask) <= 3;
    
    %Creat an FOV mask with the out pixels zeroed so we don't sample at the
    %edge (crashes 3x3 smpling)
    fov_mask = true(size(orig_im));
    fov_mask([1 end], :) = 0;
    fov_mask(:, [1 end]) = 0;
    
    %Make GT orientations
    %Create storage for the ground truth orientations
    gts_ori = zeros(size(orig_im));
    orientation_map = zeros(size(orig_im));
    
    if any(fibre_mask(:))
        
        %Extract x,y coords of vessels + vessel centres
        gts = bwmorph(fibre_mask, 'thin', inf);
        [a_y a_x] = find(fibre_mask);
        [c_y c_x] = find(gts);
        num_pts = size(c_x,1);
        
        %Loop through each skeleton point
        for jj = 1:num_pts

            %Sample local window from skeleton map
            local_win = sample_window(gts, win_size, c_y(jj), c_x(jj), 0);

            %Get all points connected to the centre
            [yi xi] = find(bwselect(local_win, (win_size+1)/2, (win_size+1)/2, 8));
            uni_x = unique(xi);
            uni_y = unique(yi);

            if length(uni_x) > length(uni_y)
                uni_y = sparse(xi, 1, yi, win_size, 1) ./ sparse(xi, 1, 1, win_size, 1);
                uni_y = full(uni_y(uni_x));
            else
                uni_x = sparse(yi, 1, xi, win_size, 1) ./ sparse(yi, 1, 1, win_size, 1);
                uni_x = full(uni_x(uni_y));
            end

            uu = mean(diff(uni_x));
            vv = -mean(diff(uni_y));
            dd = sqrt(uu^2 + vv^2);
            gts_ori(c_y(jj), c_x(jj)) = complex(uu / dd, vv / dd);
        end

        a_u = griddata(c_x, c_y, real(gts_ori(gts)),a_x, a_y, 'nearest'); %#ok
        a_v = griddata(c_x, c_y, imag(gts_ori(gts)),a_x, a_y, 'nearest'); %#ok

        dd = a_u.^2 + a_v.^2;

        orientation_map = zeros(size(orig_im));
        orientation_map(fibre_mask) = (complex(a_u, a_v).^2) ./ dd;
    end
    
    %Save stuff
    if ii <= 200
        data = 'training\';
    else
        data = 'test\';
    end
    
    base_name = fibre_list(ii).name(1:end-4);
    save([asymmetryroot 'data\fibre\' data 'images\' base_name '.mat'], 'orig_im');
    save([asymmetryroot 'data\fibre\' data 'fibre_masks\' base_name 'f_mask.mat'], 'fibre_mask');
    save([asymmetryroot 'data\fibre\' data 'red_fibre_masks\' base_name 'r_mask.mat'], 'red_fibre_mask');
    save([asymmetryroot 'data\fibre\' data 'blue_fibre_masks\' base_name 'b_mask.mat'], 'blue_fibre_mask');
    save([asymmetryroot 'data\fibre\' data 'junction_fibre_masks\' base_name 'j_mask.mat'], 'junction_mask');
    save([asymmetryroot 'data\fibre\' data 'fov_masks\' base_name 'fov_mask.mat'], 'fov_mask');
    save([asymmetryroot 'data\fibre\' data 'orientation_maps\' base_name 'ori_map.mat'], 'orientation_map');
    save([asymmetryroot 'data\fibre\' data 'fibre_masks_dilated\' base_name 'f_mask.mat'], 'fibre_mask_d3');
    
    if ii < 5
        figure; 
        a(1) = subplot(2,3,1); imgray(anno_im);
        a(2) = subplot(2,3,2); imgray(fibre_mask);
        a(3) = subplot(2,3,3); imgray(red_fibre_mask);
        a(4) = subplot(2,3,4); imgray(blue_fibre_mask);
        a(5) = subplot(2,3,5); imgray(junction_mask);
        a(6) = subplot(2,3,6); imgray(complex2rgb(orientation_map));
        linkaxes(a);
    end
end
%%
mkdir([asymmetryroot 'data\fibre\training\fibre_masks_dilated\']);
mkdir([asymmetryroot 'data\fibre\test\fibre_masks_dilated\']);
fibre_list = dir([asymmetryroot 'data\fibre\CCMdatabase_Analysed\*.bmp']);


for ii = 1:400
    
    if ii <= 200
        data = 'training\';
    else
        data = 'test\';
    end
    base_name = fibre_list(ii).name(1:end-4);
    
    load([asymmetryroot 'data\fibre\' data 'fibre_masks\' base_name 'f_mask.mat'], 'fibre_mask');
    
    fibre_mask_d3 = bwdist(fibre_mask) <= 3;
    
    if ii < 10
        figure;
        subplot(1,2,1); imgray(fibre_mask);
        subplot(1,2,2); imgray(fibre_mask_d3);
    end
    save([asymmetryroot 'data\fibre\' data 'fibre_masks_dilated\' base_name 'f_mask.mat'], 'fibre_mask_d3');
    
end
    
    
    
    
    
    


