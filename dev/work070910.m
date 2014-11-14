cd C:\isbe\dev\
load C:\isbe\dev\files\u_files.mat

er_u = model_errors3(u_files1, weights_u);
save C:\isbe\dev\scale\er_u er_u;
%%
for ii = 1:8
    for jj = 1:7
        weights = 0.5e-3*[ii jj];
        er = model_errors3(u_files1, weights);
        er_name = ['C:\isbe\dev\scale\er',...
            num2str(ii), num2str(jj)];
        save(er_name, 'er');
        clear weights er er_name;
        pack;
    end
end
%%
for ii = 1:8
    for jj = 1:7
        
        er_name = ['C:\isbe\dev\scale\er',...
            num2str(ii), num2str(jj)];
        load(er_name);
        temp.shape = zeros(101,1);
        temp.tex = zeros(101,1);
        temp.scale = zeros(101,1);
        temp.indie = zeros(101,1);
        temp.combined = zeros(101,1);

        for kk = 1:101
            temp.shape(kk) = er(kk).shape;
            temp.tex(kk) = er(kk).tex;
            temp.scale(kk) = er(kk).scale;
            temp.indie(kk) = er(kk).indie;
            temp.combined(kk) = er(kk).combined;
        end
        clear er; er = temp; clear temp;
        save(er_name, 'er'); clear er;
    end
end
%%
ers = zeros(8,7);
for ii = 1:8
    for jj = 1:7
        
        er_name = ['C:\isbe\dev\scale\er',...
            num2str(ii), num2str(jj)];
        load(er_name);
        ers(ii,jj) = mean(er.combined);
        clear er;
    end
end
%%
er_name = 'C:\isbe\dev\scale\er_r50';
er = er_r50;
nn = 50;
temp.shape = zeros(nn,1);
temp.tex = zeros(nn,1);
temp.scale = zeros(nn,1);
temp.indie = zeros(nn,1);
temp.combined = zeros(nn,1);

for kk = 1:nn
    temp.shape(kk) = er(kk).shape;
    temp.tex(kk) = er(kk).tex;
    temp.scale(kk) = er(kk).scale;
    temp.indie(kk) = er(kk).indie;
    temp.combined(kk) = er(kk).combined;
end
clear er; er_r50 = temp; clear temp;
save(er_name, 'er_r50');
