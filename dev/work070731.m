%%
N = 179;
ssv = 200;
ri = 1:179;%ceil(179*rand(1,N));
%
display('new data');
X_shape = zeros(N, 2*ssv);
X_area = zeros(N, 1);
for ii = 1:N
    load(['C:\isbe\dev\masses\', m_files(ri(ii)).name]);
    idx = round(linspace(1, size(mass.mass_outline,1), ssv+1));
    idx(end) = []; %ensures first point is not equal to the last point!
    X_shape(ii,:) = [mass.mass_outline(idx,1)', mass.mass_outline(idx,2)'];
    X_area(ii) = mass.mass_area;
    clear mass idx;
end
mean_area = mean(X_area);
%%
f(N+2) = 0;
for ii = 1:N+2; f(ii) = figure('WindowStyle', 'docked'); end

shape_one = [X_shape(1, 1:ssv)' X_shape(1, ssv+1:end)'];

figure(f(1));
plot(X_shape(1, 1:ssv), X_shape(1, ssv+1:end));
hold on;
plot(X_shape(1, 1), X_shape(1, ssv+1), 'rx');

figure(f(N+1));
plot(X_shape(1, 1:ssv), X_shape(1, ssv+1:end));
hold on;
plot(X_shape(1, 1), X_shape(1, ssv+1), 'rx');

X_pro = zeros(N, 2*ssv);
ddd = zeros(N,1);
for ii = 1:N
    figure(f(ii));
    plot(X_shape(ii, 1:ssv), X_shape(ii, ssv+1:end));
    hold on;
    plot(X_shape(ii, 1), X_shape(ii, ssv+1), 'rx');

    [dd Z t] = mb_procrustes(shape_one,...
        [X_shape(ii, 1:ssv)' X_shape(ii, ssv+1:end)']);
    X_pro(ii, :) = [Z(:,1)', Z(:,2)'];
    ddd(ii) = dd;
    %shape_scale(ii) = t.b;
    plot(Z(:,1), Z(:,2), 'g:');

    figure(f(N+1));
    plot(Z(:,1), Z(:,2), 'b:'); hold on;
    plot(Z(1,1), Z(1,2), 'rx');
    clear dd Z t;
end
mean_dd = mean(ddd);
display(['mean change = ', num2str(mean_dd)]); 
mean_shape = mean(X_pro);
%
figure(f(N+2));
plot(mean_shape(1:ssv), mean_shape(ssv+1:end), 'r'); hold on;
plot(mean_shape(1), mean_shape(ssv+1), 'bx'); hold on;

mean_shape = reshape(mean_shape, ssv, 2);
mean_shape = mean_shape - repmat(mean(mean_shape), size(mean_shape,1), 1);

scale_factor = sqrt(mean_area/polyarea(mean_shape(:,1), mean_shape(:,2)));
mean_shape = mean_shape*scale_factor;

th = atan2(mean_shape(1,1), mean_shape(1,2));
mean_shape = ([cos(th) -sin(th); sin(th) cos(th)] * mean_shape')';

plot(mean_shape(:,1), mean_shape(:,2), 'g:');
plot(mean_shape(1,1), mean_shape(1,2), 'kx');
%
%colours = 'yrcmk';
%

go_on = 1; jj = 1;
while go_on
figure(f(N+1)); hold off;

    for ii = 1:N
        [dd Z t] = mb_procrustes(mean_shape,...
            [X_pro(ii, 1:ssv)' X_pro(ii, ssv+1:end)']);
        X_pro(ii, :) = [Z(:,1)', Z(:,2)'];
        ddd(ii) = dd;
        figure(f(ii));
        plot(Z(:,1), Z(:,2), [colours(jj), ':']);
        figure(f(N+1));
        plot(Z(:,1), Z(:,2), 'b:'); hold on;
        plot(Z(1,1), Z(1,2), 'rx');
        clear dd Z t;
    end
    go_on = abs(mean_dd - mean(ddd)) > 1e-3;
    mean_dd = mean(ddd);
    display(['mean change = ', num2str(mean(ddd))]); 
    mean_shape = mean(X_pro);
    mean_shape = reshape(mean_shape, ssv, 2);
    
    figure(f(N+2));
    plot(mean_shape(:,1), mean_shape(:,2), [colours(jj), ':']);
    plot(mean_shape(1,1), mean_shape(1,2), 'kx');

    
    mean_shape = mean_shape - repmat(mean(mean_shape), size(mean_shape,1), 1);

    scale_factor = sqrt(mean_area/polyarea(mean_shape(:,1), mean_shape(:,2)));
    mean_shape = mean_shape*scale_factor;

    th = atan2(mean_shape(1,1), mean_shape(1,2));
    mean_shape = ([cos(th) -sin(th); sin(th) cos(th)] * mean_shape')';
    jj = jj+1;
end
display(['stopped after ', num2str(jj-1), ' iterations']);
%%
k_shape = length(mass_model.L_shape);
k_tex   = length(mass_model.L_tex);

W_shape = k_shape / sum(sqrt(mass_model.L_shape));
W_tex   = k_tex / sum(sqrt(mass_model.L_tex));
W_scale = 1 / sqrt(mass_model.L_scale); %length L_scale = 1
W_n = 1 / sqrt(mass_model.L_n); %length L_n = 1

mass_model.W_shape = W_shape;
mass_model.W_tex = W_tex;
mass_model.W_scale = W_scale;
mass_model.W_n = W_n;

combined_data = [W_shape*mass_model.B_shape; W_tex*mass_model.B_tex; W_scale*mass_model.B_scale]';