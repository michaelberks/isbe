%% 1) Locate a set of apexes in OCT images given manual annotations

%Get list of annotated nailfolds
nf_files = dir('C:\isbe\nailfold\images\anonymous_oct\annotations_qt\*.txt');
num_nf = length(nf_files);

mkdir C:\isbe\nailfold\playground\model_experiments\images\orig
mkdir C:\isbe\nailfold\playground\model_experiments\images\g1d
mkdir C:\isbe\nailfold\playground\model_experiments\images\g2d
    
for nn = 1:num_nf
    image_path = ['C:\isbe\nailfold\images\anonymous_oct\bmp\' nf_files(nn).name(1:end-11) '.bmp'];
    nailfold = imread(image_path);
    figure; imgray(nailfold);
    
    %Load in the annotated vessels
    v = read_vessels_from(['C:\isbe\nailfold\images\anonymous_oct\annotations_qt\' nf_files(nn).name]);
    num_vessels = length(v);
    
    %Load in the vessel apex pts
    v_apex = u_load(['C:\isbe\nailfold\images\anonymous_oct\annotations_qt\' nf_files(nn).name(1:end-4) '_v_pts.mat']);
    
    %Compute the Gaussian derivatives
    [mag_2d, ori_2d] = gaussian_2nd_derivative_line(nailfold, 4);
    [mag_1d] = gaussian_1st_derivative_gradient(nailfold, 2);
    %
    [rows cols] = size(nailfold);
    
    apexes = cell(num_vessels,1);
    
    for ii = 1:num_vessels
       
        %Discard duplicate points
        keep = [true; any(diff(v{ii}),2)];
        vessel = [v{ii}(keep,1) v{ii}(keep,2)];

        %Sample points evenly along the vessel
        dists = cumsum([0; sqrt(sum(diff(vessel).^2, 2))]);
        vessel = interp1(dists, vessel, linspace(0, dists(end), floor(dists(end))), 'spline');
    
        %Get vessel top given marked points
        if v_apex(ii,1) > v_apex(ii,2)
            vessel = flipud(vessel);
            v_apex(ii,:) = size(vessel,1) + 1 - v_apex(ii,:);
        end

        vessel_l = vessel(v_apex(ii,1):v_apex(ii,2)-1,:);
        dists = cumsum([0; sqrt(sum(diff(vessel_l).^2, 2))]);
        vessel_l = interp1(dists, vessel_l, linspace(0, dists(end), 15), 'spline');

        vessel_r = vessel(v_apex(ii,2)+1:v_apex(ii,3),:);
        dists = cumsum([0; sqrt(sum(diff(vessel_r).^2, 2))]);
        vessel_r = interp1(dists, vessel_r, linspace(0, dists(end), 15), 'spline');

        vessel_top = [vessel_l; vessel(v_apex(ii,2),:); vessel_r];
        
        plot(vessel(:,1),vessel(:,2));
        plot(vessel(v_apex(ii,2),1), vessel(v_apex(ii,2),2), 'rx');
        plot(vessel_top(:,1), vessel_top(:,2), 'm-.');
        plot(vessel_top(1,1), vessel_top(1,2), 'gx');
        plot(vessel_top(end,1), vessel_top(end,2), 'yx');

        sr = max(1, floor(min(vessel_top(:,2)) - 50));
        er = min(rows, ceil(max(vessel_top(:,2)) + 50));
        sc = max(1, floor(min(vessel_top(:,1)) - 50));
        ec = min(cols, floor(max(vessel_top(:,1)) + 50));
       
%         %Sample patch from 1st deriv.
%         mag_patch = mag_1d(sr:er, sc:ec);
% %         mag_patch = (mag_patch - min(mag_patch(:))) / (max(mag_patch(:)) - min(mag_patch(:))); %Scale patch between 0 and 1
%         save(['C:\isbe\nailfold\playground\model_experiments\images\g1d\vessel_top' zerostr(nn,3) '_' zerostr(ii,3) '.mat'], 'mag_patch');
%         imwrite(mag_patch, ['C:\isbe\nailfold\playground\model_experiments\images\g1d\vessel_top' zerostr(nn,3) '_' zerostr(ii,3) '.png']);
% 
%         %Sample patch from 1st deriv.
%         mag_patch = mag_2d(sr:er, sc:ec);
% %         mag_patch = (mag_patch - min(mag_patch(:))) / (max(mag_patch(:)) - min(mag_patch(:))); %Scale patch between 0 and 1
%         save(['C:\isbe\nailfold\playground\model_experiments\images\g2d\vessel_top' zerostr(nn,3) '_' zerostr(ii,3) '.mat'], 'mag_patch');
%         imwrite(mag_patch, ['C:\isbe\nailfold\playground\model_experiments\images\g2d\vessel_top' zerostr(nn,3) '_' zerostr(ii,3) '.png']);
% 
%         %Sample patch from 1st deriv.
%         image_patch = nailfold(sr:er, sc:ec);
%         imwrite(image_patch, ['C:\isbe\nailfold\playground\model_experiments\images\orig\vessel_top' zerostr(nn,3) '_' zerostr(ii,3) '.png']);
% 
%         %Write out a pts file we can read in to VXL
%         num_pts = size(vessel_top, 1);
%         fid1 = fopen(['C:\isbe\nailfold\playground\model_experiments\points\vessel_top' zerostr(nn,3) '_' zerostr(ii,3) '.pts'], 'wt');
%         fprintf(fid1, '%s \n', 'version: 1');
%         fprintf(fid1, '%s %d \n', 'n_points:', num_pts);
%         fprintf(fid1, '%s \n', '{'); 
%         for jj = 1:num_pts
%             fprintf(fid1,'%.2f %.2f \n', vessel_top(jj,1) - sc, vessel_top(jj,2) - sr);
%         end
%         fprintf(fid1, '%s \n', '}');
%         fprintf(fid1, 'nailfold: %s \n', image_path);
%         fprintf(fid1, '%s %d \n', 'start_row:', sr);
%         fprintf(fid1, '%s %d \n', 'start_col: ', sc);
%         fclose(fid1);
        
        %Wtie out points as a matlab file
        v_pts = [vessel_top(:,1) - sc, vessel_top(:,2) - sr];
        save(['C:\isbe\nailfold\playground\model_experiments\points_m\vessel_top' zerostr(nn,3) '_' zerostr(ii,3) '.mat'], 'v_pts', 'sr', 'sc');
        
        %Store the vessel apexes for this nailfold
        apexes{ii} = vessel_top;
    end
    save(['C:\isbe\nailfold\data\apexes\' nf_files(nn).name(1:end-4) '_apexes.mat'], 'apexes');
end
%%
%% 1) Locate a set of apexes in OCT images given manual annotations

%Get list of annotated nailfolds
nf_files = dir('C:\isbe\nailfold\images\anonymous_oct\annotations_qt\*.txt');
num_nf = length(nf_files);

mkdir C:\isbe\nailfold\playground\model_experiments\images\orig
mkdir C:\isbe\nailfold\playground\model_experiments\images\g1d
mkdir C:\isbe\nailfold\playground\model_experiments\images\g2d
vessel_num = 1;
data = 'training';

for nn = 1:num_nf
    
    if nn == 7
        data = 'test';
        vessel_num = 1;
    end
    
    image_path = ['C:\isbe\nailfold\images\anonymous_oct\bmp\' nf_files(nn).name(1:end-11) '.bmp'];
    nailfold = imread(image_path);
    figure; imgray(nailfold);
    
    %Load in the annotated vessels
    v = read_vessels_from(['C:\isbe\nailfold\images\anonymous_oct\annotations_qt\' nf_files(nn).name]);
    num_vessels = length(v);
    
    %Load in the vessel apex pts
    v_apex = u_load(['C:\isbe\nailfold\images\anonymous_oct\annotations_qt\' nf_files(nn).name(1:end-4) '_v_pts.mat']);
    
    %Compute the Gaussian derivatives
    [mag_2d, ori_2d] = gaussian_2nd_derivative_line(nailfold, 4);
    [mag_1d] = gaussian_1st_derivative_gradient(nailfold, 2);
    %
    [rows cols] = size(nailfold);
    
    apexes = cell(num_vessels,1);
    
    for ii = 1:num_vessels
       
        %Discard duplicate points
        keep = [true; any(diff(v{ii}),2)];
        vessel = [v{ii}(keep,1) v{ii}(keep,2)];

        %Sample points evenly along the vessel
        dists = cumsum([0; sqrt(sum(diff(vessel).^2, 2))]);
        vessel = interp1(dists, vessel, linspace(0, dists(end), floor(dists(end))), 'spline');
    
        %Get vessel top given marked points
        if v_apex(ii,1) > v_apex(ii,2)
            vessel = flipud(vessel);
            v_apex(ii,:) = size(vessel,1) + 1 - v_apex(ii,:);
        end

        vessel_l = vessel(v_apex(ii,1):v_apex(ii,2)-1,:);
        dists = cumsum([0; sqrt(sum(diff(vessel_l).^2, 2))]);
        vessel_l = interp1(dists, vessel_l, linspace(0, dists(end), 15), 'spline');

        vessel_r = vessel(v_apex(ii,2)+1:v_apex(ii,3),:);
        dists = cumsum([0; sqrt(sum(diff(vessel_r).^2, 2))]);
        vessel_r = interp1(dists, vessel_r, linspace(0, dists(end), 15), 'spline');

        vessel_top = [vessel_l; vessel(v_apex(ii,2),:); vessel_r];
        
        plot(vessel(:,1),vessel(:,2));
        plot(vessel(v_apex(ii,2),1), vessel(v_apex(ii,2),2), 'rx');
        plot(vessel_top(:,1), vessel_top(:,2), 'm-.');
        plot(vessel_top(1,1), vessel_top(1,2), 'gx');
        plot(vessel_top(end,1), vessel_top(end,2), 'yx');

        sr = max(1, floor(min(vessel_top(:,2)) - 100));
        er = min(rows, ceil(max(vessel_top(:,2)) + 100));
        sc = max(1, floor(min(vessel_top(:,1)) - 100));
        ec = min(cols, floor(max(vessel_top(:,1)) + 100));

        %Sample patch from 1st deriv.
        image_patch = nailfold(sr:er, sc:ec);
        imwrite(image_patch, ['C:\isbe\nailfold\data\' data '\patches\apex' zerostr(vessel_num, 4) '.png']);

        %Write out a pts file we can read in to VXL
        num_pts = size(vessel_top, 1);
        fid1 = fopen(['C:\isbe\nailfold\data\' data '\patch_apexes\apex' zerostr(vessel_num, 4) '.pts'], 'wt');
        fprintf(fid1, '%s \n', 'version: 1');
        fprintf(fid1, '%s %d \n', 'n_points:', num_pts);
        fprintf(fid1, '%s \n', '{'); 
        for jj = 1:num_pts
            fprintf(fid1,'%.2f %.2f \n', vessel_top(jj,1) - sc, vessel_top(jj,2) - sr);
        end
        fprintf(fid1, '%s \n', '}');
        fprintf(fid1, 'nailfold: %s \n', image_path);
        fprintf(fid1, '%s %d \n', 'start_row:', sr);
        fprintf(fid1, '%s %d \n', 'start_col: ', sc);
        fclose(fid1);
        
        %Wtie out points as a matlab file
        v_pts = [vessel_top(:,1) - sc, vessel_top(:,2) - sr];
        save(['C:\isbe\nailfold\data\' data '\patch_apexes\apex' zerostr(vessel_num, 4) '.mat'], 'v_pts', 'sr', 'sc');
        
        %Icrement the vessel count
        vessel_num = vessel_num + 1;
    end
    
end
%%
%Make masks for the images
nf_files = dir('C:\isbe\nailfold\images\anonymous_oct\annotations_qt\*.txt');
num_nf = length(nf_files);

vessel_num = 1;
data = 'training';

mkdir('C:\isbe\nailfold\images\anonymous_oct\masks\');
for nn = 1:num_nf
    
    if nn == 7
        data = 'test';
        vessel_num = 1;
    end
    
    image_path = ['C:\isbe\nailfold\images\anonymous_oct\bmp\' nf_files(nn).name(1:end-11) '.bmp'];
    mask_path = ['C:\isbe\nailfold\data\' data '\images\nailfold' zerostr(nn,3) '_mask.bmp'];
    
    nailfold = imread(image_path);
    
    nailfold_mask = nailfold > 240;
    edges = false(size(nailfold));
    edges([1 end], :) = 1;
    edges(:, [1 end]) = 1;
    [edge_r edge_c] = find(edges);
    
    nailfold_mask = ~bwselect(nailfold_mask, edge_c, edge_r);
    nailfold_mask = imerode(nailfold_mask, strel('disk', 20));
    nailfold_mask = uint8(nailfold_mask*255);
    imwrite(nailfold_mask, mask_path);
    save(['C:\isbe\nailfold\images\anonymous_oct\masks\' nf_files(nn).name(1:end-11) '_mask.mat'], 'nailfold_mask');
    figure;
    subplot(2,1,1); imgray(nailfold);
    subplot(2,1,2); imgray(nailfold_mask);
end
%%
nf_files = dir('C:\isbe\nailfold\images\anonymous_oct\annotations_qt\*.txt');
num_nf = length(nf_files);

mkdir C:\isbe\nailfold\data\training\apexes\
mkdir C:\isbe\nailfold\data\test\apexes\
data = 'training';

for nn = 1:num_nf
    
    if nn == 7
        data = 'test';
        vessel_num = 1;
    end
    
    apexes_path = ['C:\isbe\nailfold\data\apexes\' nf_files(nn).name(1:end-4) '_apexes.mat'];
    new_path = ['C:\isbe\nailfold\data\' data '\apexes\nailfold' zerostr(nn,3) '_apexes.mat'];
    copyfile(apexes_path, new_path);
end