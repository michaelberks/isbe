load C:\isbe\dev\files\small_files.mat

N = length(small_files);
size_shape_vec = 500;
shapes_ori = zeros(N, 2*size_shape_vec);
shapes = zeros(2, size_shape_vec, N);
mass_areas = zeros(N,1);

for ii = 1:N
    temp = load(['C:\isbe\dev\masses\', small_files(ii).name]);
    mass = temp.mass; clear temp;
    
    shape_vec = mass.mass_outline;
    idx = round(linspace(1, length(shape_vec(:,1)), size_shape_vec+1));
    idx = idx(1:size_shape_vec); %ensures first point is not equal to the last point!
    shapes_ori(ii,:) = [shape_vec(idx,1)', shape_vec(idx,2)'];    
    mass_areas(ii) = mass.mass_area;
    clear mass;
end

[mean_X, P_x, B_x, L_x,...
    mean_scale, P_scale, B_scale, L_scale mean_target, shape_scale, shapes_pro]...
    = shape_model(shapes_ori, mean(mass_areas), 0.98);
clear mean_X P_x B_x L_x mean_scale P_scale B_scale L_scale mean_target shape_scale;
%%
for ii = 1:N
    shapes(:,:,ii) = reshape(shapes_ori(ii,:),[],2)';
end
%%
N = 25;%length(small_files);
size_shape_vec = 500;
%X_shape = zeros(N, 2*size_shape_vec);
shapes = zeros(2, size_shape_vec, N);
mass_areas = zeros(N,1);

for ii = 1:N
    temp = load(['C:\isbe\dev\masses\', small_files(ii).name]);
    mass = temp.mass; clear temp;
    
    shape_vec = mass.mass_outline;
    idx = round(linspace(1, length(shape_vec(:,1)), size_shape_vec+1));
    idx(end) = []; %ensures first point is not equal to the last point!
    shapes(:,:,ii) = [shape_vec(idx,1)'; shape_vec(idx,2)'];    
    
    clear mass;
end
%%
profile on;
[reparameterisedPoints, params, objValRecord, nFevalsRecord, timeTakenRecord, outputFileName] = amb_automodelbuild (shapes, ...
   'circle', 'nIterations', 100, 'saveFrequency', 10, 'Quiet', 0, 'optimisePose', 1, 'optimiseOrigin', 0);
profile viewer;
% the results are saved into a file, which can be summarised into a HTML
% file using the following command:
%amb_generate_report (outputFileName);
%%
for ii = 1:10
    figure('WindowStyle', 'docked'); hold on;
    plot(reparameterisedPoints(1,:,ii), reparameterisedPoints(2,:,ii), 'r-x');
    plot(shapes(1,:,ii), shapes(2,:,ii), 'b-.');
end
%%
figure('WindowStyle', 'docked'); hold on;
for ii = 1:25
    plot(reparameterisedPoints(1,:,ii), reparameterisedPoints(2,:,ii), 'r-x');
end
%%
figure('WindowStyle', 'docked'); hold on;
for ii = 1:20
    plot(shapes_m(1,:,ii), shapes_m(2,:,ii), 'b-.');
end
%%
for ii = 1:50
    file_out = ['C:\isbe\dev\models\loo_r50\model', zerostr(ii, 3)];
    one_out_files = r_files50;
    one_out_files(ii) = [];

    mass_model = generate_mass_AM(...
        one_out_files, file_out, 500, 'C:\isbe\dev\masses\', 20000);
end