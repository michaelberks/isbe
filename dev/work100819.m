%%
%%
slice = cell(13,6);
data = cell(12,7);
letters1 = 'MABCDEFGHIJKL';
letters2 = 'MJIHGFEDCBALK';
for ii = 1:13
    for jj = 1:6
        slice{ii,jj} = [letters2(ii) num2str(jj)];
    end
end
%
data(:,1) = [slice(1,1:6)'; slice(1,1:6)'];
for ii = 1:12;
    data(ii,2:7) = circshift(slice(1+ii,:), [0 1-ii]);
end
for ii = 1:6;
    data(:,1+ii) = circshift(data(:,1+ii), [ii-1 0]);
end
data(7:12,:) = conjstr(data(7:12,:))
%%
slice = cell(3,13,6,2);
data = cell(3,12,7,2);
letters2 = 'MABCDEFGHIJKL';
%letters2 = 'MJIHGFEDCBALK';
for kk = 1:3
    for ii = 1:13
        for jj = 1:6
            for ll = 1:2
                slice{kk,ii,jj,ll} = [letters2(ii) num2str(jj) '(' num2str(kk) ',' num2str(ll) ')'];
            end
        end
    end
end
%
data(:,1:6,1,:) = permute(slice(:,1,1:6,:), [1 3 2 4]);
data(:,7:12,1,:) = permute(slice(:,1,1:6,:), [1 3 2 4]);
%
for ii = 1:12;
    data(:,ii,2:7,:) = circshift(slice(:,1+ii,:,:), [0 0 1-ii 0]);
end
%
for ii = 2:6;
    data(:,:,1+ii,:) = circshift(data(:,:,1+ii,:), [0 ii-1 0 0]);
end
data(:,7:12,:,:) = conjstr(data(:,7:12,:,:));

display(squeeze(data(1,:,:,:)))
display(squeeze(data(2,:,:,:)))
display(squeeze(data(3,:,:,:)))
%%
%--------------------------------------------------------------------------
bg = u_load('C:\isbe\asymmetry_project\data\synthetic_backgrounds\smooth128\bg001.mat');
[cls label] = create_gauss_bar(3, 8, 23, 128, 128, 64, 64);
cls = bg+cls;
figure; imagesc(cls); axis image; hold on;

n = 12;
r = [60 60 80 80]';
c = [60 80 80 60]';
rr = [r bsxfun(@plus, r, sin(2*(3:-1:-8)*pi/n))];
cc = [c bsxfun(@plus, c, cos(2*(3:-1:-8)*pi/n))];

colors = lines(6);
for ii = 1:6
    plot(cc(:,ii+[1 7]), rr(:,ii+[1 7]), 'markeredgecolor', colors(ii,:), 'marker', 'o', 'linestyle', 'none');
end
letters = 'MABCDEFGHIJKLM';
for ii = 1:13
    text(cc(:,ii), rr(:,ii), letters(ii));
end

dt = dtwavexfm2b(cls, 2);
slice = dt_to_pixel_subset(dt, rr, cc);
data = zeros(4,12,7,2);

data(:,1:6,1,:) = permute(slice(:,1,1:6,:), [1 3 2 4]);
data(:,7:12,1,:) = permute(slice(:,1,1:6,:), [1 3 2 4]);
%
for ii = 1:12;
    data(:,ii,2:7,:) = circshift(slice(:,1+ii,:,:), [0 0 1-ii 0]);
end
%
for ii = 2:6;
    data(:,:,1+ii,:) = circshift(data(:,:,1+ii,:), [0 ii-1 0 0]);
end
data(:,7:12,:,:) = conj(data(:,7:12,:,:));

display(squeeze(data(1,:,:,:)))
display(squeeze(data(2,:,:,:)))
display(squeeze(data(3,:,:,:)))
%%
%**************************************************************************
%**************************************************************************
slice = cell(3,7,6,2);
data = cell(3,6,7,2);
letters2 = 'MABCDEF';
for kk = 1:3
    for ii = 1:7
        for jj = 1:6
            for ll = 1:2
                slice{kk,ii,jj,ll} = [letters2(ii) num2str(jj) '(' num2str(kk) ',' num2str(ll) ')'];
            end
        end
    end
end
%

%
for ii = 1:6;
    data(:,ii,2:7,:) = circshift(slice(:,1+ii,:,:), [0 0 2*(1-ii) 0]);
end
data = permute(data, [1 3 2 4]);
data(:,1,:,:) = slice(:,1,:,:); 

display(squeeze(data(1,:,:,:)))
display(squeeze(data(2,:,:,:)))
display(squeeze(data(3,:,:,:)))
%%
%--------------------------------------------------------------------------
bg = u_load('C:\isbe\asymmetry_project\data\synthetic_backgrounds\smooth128\bg001.mat');
[cls label] = create_gauss_bar(3, 8, 23, 128, 128, 64, 64);
cls = bg+cls;
figure; imagesc(cls); axis image; hold on;

n = 6;
r = [60 60 80 80]';
c = [60 80 80 60]';
rr = [r bsxfun(@plus, r, sin(2*(0:-1:-5)*pi/n))];
cc = [c bsxfun(@plus, c, cos(2*(0:-1:-5)*pi/n))];

colors = lines(6);
for ii = 1:6
    plot(cc(:,ii+1), rr(:,ii+1), 'markeredgecolor', colors(ii,:), 'marker', 'o', 'linestyle', 'none');
end
letters = 'MABCDEF';
for ii = 1:7
    text(cc(:,ii), rr(:,ii), letters(ii));
end
%%
dt = dtwavexfm2b(cls, 2);
slice = dt_to_pixel_subset(dt, rr, cc);
data = zeros(4,6,7,2);

data(:,:,1,:) = permute(slice(:,1,:,:), [1 3 2 4]);
%
for ii = 1:6;
    data(:,ii,2:7,:) = circshift(slice(:,1+ii,:,:), [0 0 2*(1-ii) 0]);
end

display(squeeze(data(1,:,:,:)))
display(squeeze(data(2,:,:,:)))
display(squeeze(data(3,:,:,:)))