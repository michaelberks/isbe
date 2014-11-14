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
    
    for ii = 1:length(v)
       
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
    end
end
%%

    
%%
    
%% 3) Some testing stuff...
test_names = [];
for ii = 7:12
    
    for jj = 1:30
        if exist(['C:\isbe\nailfold\playground\model_experiments\points\vessel_top' sprintf('%03d', ii) '_'  sprintf('%03d', jj) '.pts'], 'file');
            test_names{end+1,1} = [sprintf('%03d', ii) '_'  sprintf('%03d', jj)]; %#ok
        else
            break;
        end
    end
end
train_names = [];
for ii = 1:6
    
    for jj = 1:30
        if exist(['C:\isbe\nailfold\playground\model_experiments\points\vessel_top' sprintf('%03d', ii) '_'  sprintf('%03d', jj) '.pts'], 'file');
            train_names{end+1,1} = [sprintf('%03d', ii) '_'  sprintf('%03d', jj)]; %#ok
        else
            break;
        end
    end
end
save('C:\isbe\nailfold\playground\model_experiments\vessel_names.mat', 'test_names', 'train_names');
%
feature = 'orig';

num_rows = 3;
num_cols = 4;
jj = 1;
for ff = 1:10
    figure;
    for row = 1:num_rows, 
        for col = 1:num_cols        

            f1 = fopen(['C:\isbe\nailfold\playground\model_experiments\points\vessel_top' test_names{jj} '.pts']);
            textscan(f1, '%[^{]');
            fgetl(f1);
            vessel_str = textscan(f1,'%[^}]');
            fclose(f1);
            orig_pts = str2num(vessel_str{1}{1});

            f1 = fopen(['C:\isbe\nailfold\playground\model_experiments\out_points\' feature '\vessel_top' test_names{jj} '.pts']);
            textscan(f1, '%[^{]');
            fgetl(f1);
            vessel_str = textscan(f1,'%[^}]');
            fclose(f1);
            test_pts = str2num(vessel_str{1}{1});
            
            f1 = fopen(['C:\isbe\nailfold\playground\model_experiments\out_points\g2d\vessel_top' test_names{jj} '.pts']);
            textscan(f1, '%[^{]');
            fgetl(f1);
            vessel_str = textscan(f1,'%[^}]');
            fclose(f1);
            test_pts2 = str2num(vessel_str{1}{1});

            im = imread(['C:\isbe\nailfold\playground\model_experiments\images\' feature '\vessel_top' test_names{jj} '.png']);

            axes('units', 'normalized', 'position', [(col-1)/num_cols (row-1)/num_rows 1/num_cols 1/num_rows]);
            imagesc(im); axis image off; colormap(gray(256));
            hold on;
            plot(orig_pts(:,1), orig_pts(:,2), 'rx');
            plot(test_pts2(:,1), test_pts2(:,2), 'cx');
            plot(test_pts(:,1), test_pts(:,2), 'gx');
            jj = jj + 1;
        end
    end
end

