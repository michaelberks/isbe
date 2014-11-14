%% 1) Locate a set of apexes in 1 image given manual annotations
% We'll use MB's markup of the KK demo image (because it just happens to be
% an image I've marked up lots of vessels on)

%load in the nailfold - note we need to resize this by a factor 3/4 to
%match the other nailfolds in our set
nailfold = imresize(imread([nailfoldroot 'images\kk_demo\nailfolds\n3_mb.bmp']), 0.75, 'bilinear');

%Load in the annotated vessels
v = u_load([nailfoldroot 'images\kK_demo\nailfolds\n3_mb_vessels_good.mat']);
vessels = [];
for ii = 1:length(v)
    if size(v{ii},1) > 1;       
        %Discard duplicate points
        keep = [true; any(diff(v{ii}),2)];
        vessel = 0.75*[v{ii}(keep,1) v{ii}(keep,2)];

        %Sample points evenly along the vessel
        dists = cumsum([0; sqrt(sum(diff(vessel).^2, 2))]);
        vessels{end+1,1} = interp1(dists, vessel, linspace(0, dists(end), floor(dists(end))), 'spline'); %#ok
    end
end
num_vessels = length(vessels);
%%
%Now choose each vessel apex as the point with highest y-val (I know this
%may not strictly be accurate given local deviations, rotations etc - but
%it appears to work more robustly than more complex methods such as
%searching for points of maximum curvature etc)
figure; hold on; axis equal ij;
vessel_apex = zeros(num_vessels, 2);
vessel_tops = cell(num_vessels, 1);
for ii = 1:num_vessels
    [y_min min_idx] = min(vessels{ii}(:,2));
    vessel_apex(ii,:) = [vessels{ii}(min_idx,1) y_min];
    
    vessel_l = vessels{ii}(min_idx:-1:1,:);
    dists = cumsum([0; sqrt(sum(diff(vessel_l).^2, 2))]);
    vessel_l = interp1(dists, vessel_l, 2:2:30, 'spline');
    
    vessel_r = vessels{ii}(min_idx:end,:);
    dists = cumsum([0; sqrt(sum(diff(vessel_r).^2, 2))]);
    vessel_r = interp1(dists, vessel_r, 2:2:30, 'spline');
    
    vessel_tops{ii} = [flipud(vessel_l); vessels{ii}(min_idx,:); vessel_r];
    if vessel_tops{ii}(1,1) > vessel_tops{ii}(end,1)
        vessel_tops{ii} = flipud(vessel_tops{ii});
    end
    
    
    plot(vessels{ii}(:,1), vessels{ii}(:,2)); 
    plot(vessel_apex(ii,1), vessel_apex(ii,2), 'rx');
    %plot(vessel_l(:,1), vessel_l(:,2), 'g.');
    %plot(vessel_r(:,1), vessel_r(:,2), 'y.');
    plot(vessel_tops{ii}(:,1), vessel_tops{ii}(:,2), 'm-.');
    plot(vessel_tops{ii}(1,1), vessel_tops{ii}(1,2), 'gx');
    plot(vessel_tops{ii}(end,1), vessel_tops{ii}(end,2), 'yx');
end

%% 2) Now compute Gaussian 1st + 2nd derivatives of the nailfold and sample
% a patch of each about each apex
mkdir C:\isbe\nailfold\playground\model_experiments\images\orig
mkdir C:\isbe\nailfold\playground\model_experiments\images\g1d
mkdir C:\isbe\nailfold\playground\model_experiments\images\g2d
%Compute the Gaussian derivatives
[mag_2d, ori_2d] = gaussian_clover_line(nailfold, 4);
[mag_1d] = gaussian_1st_derivative_gradient2(nailfold, 2);
%
for ii = 1:num_vessels
    sr = floor(min(vessel_tops{ii}(:,2)) - 50);
    er = ceil(max(vessel_tops{ii}(:,2)) + 50);
    sc = floor(min(vessel_tops{ii}(:,1)) - 50);
    ec = floor(max(vessel_tops{ii}(:,1)) + 50);
    
    %Sample patch from 1st deriv.
    mag_patch = mag_1d(sr:er, sc:ec);
    mag_patch = (mag_patch - min(mag_patch(:))) / (max(mag_patch(:)) - min(mag_patch(:))); %Scale patch between 0 and 1
    imwrite(mag_patch, ['C:\isbe\nailfold\playground\model_experiments\images\g1d\vessel_top' zerostr(ii,3) '.png']);
    
    %Sample patch from 1st deriv.
    mag_patch = mag_2d(sr:er, sc:ec);
    mag_patch = (mag_patch - min(mag_patch(:))) / (max(mag_patch(:)) - min(mag_patch(:))); %Scale patch between 0 and 1
    imwrite(mag_patch, ['C:\isbe\nailfold\playground\model_experiments\images\g2d\vessel_top' zerostr(ii,3) '.png']);
    
    %Sample patch from 1st deriv.
    image_patch = nailfold(sr:er, sc:ec);
    imwrite(image_patch, ['C:\isbe\nailfold\playground\model_experiments\images\orig\vessel_top' zerostr(ii,3) '.png']);
    
    num_pts = size(vessel_tops{ii}, 1);
    fid1 = fopen(['C:\isbe\nailfold\playground\model_experiments\points\vessel_top' zerostr(ii,3) '.pts'], 'wt');
    fprintf(fid1, '%s \n', 'version: 1');
    fprintf(fid1, '%s %d \n', 'n_points:', num_pts);
    fprintf(fid1, '%s \n', '{'); 
    for jj = 1:num_pts
        fprintf(fid1,'%.2f %.2f \n', vessel_tops{ii}(jj,1) - sc, vessel_tops{ii}(jj,2) - sr);
    end
    fprintf(fid1, '%s \n', '}');    
    fclose(fid1);
end
%% 3) Some testing stuff...

for ii = 25:34
    f1 = fopen(['C:\isbe\nailfold\playground\model_experiments\points\vessel_top' sprintf('%03d', ii) '.pts']);
    textscan(f1, '%[^{]');
    fgetl(f1);
    vessel_str = textscan(f1,'%[^}]');
    fclose(f1);
    orig_pts = str2num(vessel_str{1}{1});
    
    f1 = fopen(['C:\isbe\nailfold\playground\model_experiments\out_points\vessel_top' sprintf('%03d', ii) '.pts']);
    textscan(f1, '%[^{]');
    fgetl(f1);
    vessel_str = textscan(f1,'%[^}]');
    fclose(f1);
    test_pts = str2num(vessel_str{1}{1});
    
    im = imread(['C:\isbe\nailfold\playground\model_experiments\images\g2d\vessel_top' sprintf('%03d', ii) '.png']);
    
    figure; imagesc(im); axis image; colormap(gray(256));
    hold on;
    plot(orig_pts(:,1), orig_pts(:,2), 'rx');
    plot(test_pts(:,1), test_pts(:,2), 'gx');
end

