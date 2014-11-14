figure('windowstyle', 'docked'); hold on;

for ii = 1:50
    plot(idx_small(ii), er_small.combined(ii), 'rx');
    plot(idx_r50(ii), er_r50.combined(ii), 'b+');
end
for ii = 1:51
    plot(idx_large(ii), er_large.combined(ii), 'gx');
    plot(idx_r51(ii), er_r51.combined(ii), 'c+');
end
for ii = 1:101
    plot(idx1(ii), er_u.combined(ii), 'mo');
end
%%
ers_combined = nan;
ers_combined = ers_combined(ones(179,3));
ers_shape = ers_combined;
ers_tex = ers_combined;
ers_scale = ers_combined;

for ii = 1:50
    ers_combined(idx_small(ii),1) = er_small.combined(ii);
    ers_shape(idx_small(ii),1) = er_small.shape(ii);
    ers_tex(idx_small(ii),1) = er_small.tex(ii);
    ers_scale(idx_small(ii),1) = er_small.scale(ii);
end
for ii = 1:51
    ers_combined(idx_large(ii), 2) = er_large.combined(ii);
    ers_shape(idx_large(ii), 2) = er_large.shape(ii);
    ers_tex(idx_large(ii), 2) = er_large.tex(ii);
    ers_scale(idx_large(ii), 2) = er_large.scale(ii);
end
for ii = 1:101
    ers_combined(idx_u1(ii), 3) = er_u.combined(ii);
    ers_shape(idx_u1(ii), 3) = er_u.shape(ii);
    ers_tex(idx_u1(ii), 3) = er_u.tex(ii);
    ers_scale(idx_u1(ii), 3) = er_u.scale(ii);
end
ers_combined = sortrows(ers_combined, 3);
ers_combined(102:179,:) = [];
ers_shape = sortrows(ers_shape, 3);
ers_shape(102:179,:) = [];
ers_tex = sortrows(ers_tex, 3);
ers_tex(102:179,:) = [];
ers_scale = sortrows(ers_scale, 3);
ers_scale(102:179,:) = [];

figure('windowstyle', 'docked'); hold on; plot(ers_combined, 'x');
figure('windowstyle', 'docked'); hold on; plot(ers_shape, 'x');
figure('windowstyle', 'docked'); hold on; plot(ers_tex, 'x');
figure('windowstyle', 'docked'); hold on; plot(ers_scale, 'x');
%%
for ii = 1:51
    temp.shape(ii) = er_large(ii).shape;
    temp.tex(ii) = er_large(ii).tex;
    temp.scale(ii) = er_large(ii).scale;
    temp.indie(ii) = er_large(ii).indie;
    temp.combined(ii) = er_large(ii).combined;
end
er_large = temp; clear temp;

for ii = 1:50
    temp.shape(ii) = er_small(ii).shape;
    temp.tex(ii) = er_small(ii).tex;
    temp.scale(ii) = er_small(ii).scale;
    temp.indie(ii) = er_small(ii).indie;
    temp.combined(ii) = er_small(ii).combined;
end
er_small = temp; clear temp;
%%
for ii = 1:101
    temp.shape(ii) = er_u(ii).shape;
    temp.tex(ii) = er_u(ii).tex;
    temp.scale(ii) = er_u(ii).scale;
    temp.indie(ii) = er_u(ii).indie;
    temp.combined(ii) = er_u(ii).combined;
end
er_u = temp; clear temp;
%%
figure('windowstyle', 'docked'); hold on;

er_temp2 = nan;
er_temp2 = er_temp(ones(179,3));
for ii = 1:50
    er_temp2(idx_r50(ii),1) = er_r50.combined(ii);
    er_temp2(idx_r50(ii),2) = er_u.combined(ii);
end

plot(er_temp, 'x');
%%
figure('windowstyle', 'docked'); hold on;

load C:\isbe\dev\weights\er3_c.mat
load C:\isbe\dev\weights\er3_o.mat
load C:\isbe\dev\weights\er3_s.mat
load C:\isbe\dev\weights\er3_z.mat

er_temp = [er3_o.combined er3_s.combined er3_c.combined er3_z.combined];
er_temp = sortrows(er_temp, 1);

plot(er_temp, 'x');
%%
figure('windowstyle', 'docked');
plot(er3_c.combined, er3_s.combined, 'rx');
figure('windowstyle', 'docked');
plot(er3_c.combined, er3_z.combined, 'rx');
figure('windowstyle', 'docked');
plot(er3_s.combined, er3_z.combined, 'rx');
figure('windowstyle', 'docked');
plot(er3_o.combined, er3_s.combined, 'rx');
    
    
    