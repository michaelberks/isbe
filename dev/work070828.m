%%
er_small = model_errors2('models\small_model', small_files, 'masses\'); save scale\er_small er_small;
er_area = model_errors2('models\area_model', area_files, 'masses\'); save scale\er_area er_area;

%%
idx_CC = [];
idx_ML = [];
for ii = 1:179
    if strcmp(m_files(ii).name(9:10), 'CC')
        idx_CC(end+1) = ii;
    else
        idx_ML(end+1) = ii;
    end
end
%%        
for jj = 1:141
    unique_files(jj).idx = [];
    for ii = 1:179
        if strcmp(m_files(ii).name(5:7), zerostr(jj,3)) 
            unique_files(jj).idx(end +1) = ii;
        end
    end 
end

%%
for ii = 1:179;
    masses(ii).name = m_files(ii).name([5:8, 11]);
    display(masses(ii).name);
end
%%
clear idx_u
idx_u(1) = 1;
for ii = 2:179
    flags = zeros(ii-1,1);
    for jj = 1:ii-1
        flags(jj) = strcmp(masses(ii).name, masses(jj).name);
    end
    if ~sum(flags)
        idx_u(end+1) = ii;
        display(masses(ii).name);
    end
end
%%
clear unique_masses
for ii = 1:101
    unique_masses(ii).idx = [];
end
u_masses = masses(idx_u);
for ii = 1:179
    for jj = 1:101
        if strcmp(masses(ii).name, u_masses(jj).name)
            unique_masses(jj).idx(end+1) = ii;
        end
    end
end
%%
for ii = 1:101
    mstr = [];
    for jj = 1:length(unique_masses(ii).idx)
        mstr = [mstr, m_files(unique_masses(ii).idx(jj)).name, '; '];
    end
    display(mstr)
end
%%
frog = er3_o;

temp.shape = zeros(179,1);
temp.tex = zeros(179,1);
temp.scale = zeros(179,1);
temp.indie = zeros(179,1);
temp.combined = zeros(179,1);
for ii = 1:179
    temp.shape(ii) = frog(ii).shape;
    temp.tex(ii) = frog(ii).tex;
    temp.scale(ii) = frog(ii).scale;
    temp.indie(ii) = frog(ii).indie;
    temp.combined(ii) = frog(ii).combined;
end
er3_o = temp; clear temp frog;
save C:\isbe\dev\weights\er3_o er3_o
%%

figure('WindowStyle', 'docked');
plot(er3_s.combined, er_s.combined, 'rx');
figure('WindowStyle', 'docked');
plot(er3_c.combined, er_c.combined, 'rx');
figure('WindowStyle', 'docked');
plot(er3_z.combined, er_z.combined, 'rx');
figure('WindowStyle', 'docked');
plot(er3_o.combined, er_o.combined, 'rx');
figure('WindowStyle', 'docked');
plot(er3_o2.combined, er_o2.combined, 'rx');
figure('WindowStyle', 'docked');
plot(er3_o3.combined, er_o3.combined, 'rx');

%%
for ii = 1:179
    load(['C:\isbe\dev\masses\left_out_models\model', zerostr(ii, 3)]);
    load(['C:\isbe\dev\masses\', m_files(ii).name]);
    k_shape(ii)     = size(mass_model.P_shape, 2);
    k_tex(ii)       = size(mass_model.P_tex, 2);
    k_c1(ii)        = size(mass_model.P_c, 2);
    w1(ii, 1)       = mass_model.W_shape;
    w1(ii, 2)       = mass_model.W_tex;
    w1(ii, 3)       = mass_model.W_scale;
    mass_model = combine_model(mass_model, [weights_s, 1]);
    k_c2(ii)        = size(mass_model.P_c, 2);
    w2(ii, 1)       = mass_model.W_shape;
    w2(ii, 2)       = mass_model.W_tex;
    w2(ii, 3)       = mass_model.W_scale;
    clear mass mass_model;
end