%       Testing flow
%--------------------------------------------------------------------------
frames_dir = 'N:\Nailfold Capillaroscopy\wellcome\flow_data\';
flow_dir = 'C:\isbe\nailfold\data\wellcome_study\flow_results\';
frames_list = dir([frames_dir '*v01*.mat']);
num_vessels = length(frames_list);

mean_flow_erros = zeros(num_vessels,3);
mean_flow_velocity_erros = zeros(num_vessels,3);
mean_flow_direction_erros = zeros(num_vessels,3);
var_flow_velocity_erros = zeros(num_vessels,3);
var_flow_direction_erros = zeros(num_vessels,3);

results_path = 'C:\isbe\nailfold\data\wellcome_study\results\flow_validity_exp.mat';
%%
for i_ve = 3959:num_vessels
    display(['Processing vessel ' num2str(i_ve)]);
    
    if ~exist([flow_dir frames_list(i_ve).name], 'file')
        continue;
    end
    
    f = load([frames_dir frames_list(i_ve).name], 'frames_i');
    f1 = u_load([flow_dir frames_list(i_ve).name]);
%     [f1.flowPyramidEst, f1.flowConfidence] =...
%         estimate_flow_multilevel(f.frames_i(:,:,1:120), [], [], 1:3);
    [f1a.flowPyramidEst, f1a.flowConfidence] =...
        estimate_flow_multilevel(f.frames_i(:,:,1:2:end), [], [], 1:3);
    [f1b.flowPyramidEst, f1b.flowConfidence] =...
        estimate_flow_multilevel(f.frames_i(:,:,2:2:end), [], [], 1:3);
    
    % f1 vs f1a
    mean_flow_erros(i_ve,1) = mean(f1.flowPyramidEst{1}(:)-0.5*f1a.flowPyramidEst{1}(:));
    velocity_error = abs(f1.flowPyramidEst{1}(:))-abs(0.5*f1a.flowPyramidEst{1}(:));
    direction_error = angle(f1.flowPyramidEst{1}(:).*conj(0.5*f1a.flowPyramidEst{1}(:)));
    mean_flow_velocity_erros(i_ve,1) = mean(velocity_error);
    mean_flow_direction_erros(i_ve,1) = mean(direction_error);
    var_flow_velocity_erros(i_ve,1) = var(velocity_error);
    var_flow_direction_erros(i_ve,1) = var(direction_error);
    
    % f1 vs f1b
    mean_flow_erros(i_ve,2) = mean(f1.flowPyramidEst{1}(:)-0.5*f1b.flowPyramidEst{1}(:));
    velocity_error = abs(f1.flowPyramidEst{1}(:))-abs(0.5*f1b.flowPyramidEst{1}(:));
    direction_error = angle(f1.flowPyramidEst{1}(:).*conj(0.5*f1b.flowPyramidEst{1}(:)));
    mean_flow_velocity_erros(i_ve,2) = mean(velocity_error);
    mean_flow_direction_erros(i_ve,2) = mean(direction_error);
    var_flow_velocity_erros(i_ve,2) = var(velocity_error);
    var_flow_direction_erros(i_ve,2) = var(direction_error);
    
    % f1a vs f1b
    mean_flow_erros(i_ve,3) = 0.5*mean(f1a.flowPyramidEst{1}(:)-f1b.flowPyramidEst{1}(:));
    velocity_error = 0.5*abs(f1a.flowPyramidEst{1}(:))-abs(f1b.flowPyramidEst{1}(:));
    direction_error = angle(f1a.flowPyramidEst{1}(:).*conj(f1b.flowPyramidEst{1}(:)));
    mean_flow_velocity_erros(i_ve,3) = mean(velocity_error);
    mean_flow_direction_erros(i_ve,3) = mean(direction_error);
    var_flow_velocity_erros(i_ve,3) = var(velocity_error);
    var_flow_direction_erros(i_ve,3) = var(direction_error);
    
    if ~rem(i_ve, 100)
        save(results_path, ...
            'mean_flow_erros',...
            'mean_flow_velocity_erros',...
            'mean_flow_direction_erros',...
            'var_flow_velocity_erros',...
            'var_flow_direction_erros');
    end
    
    if i_ve <= 10
        figure; 
        subplot(2,3,1); show_flow_as('rgb', f1.flowPyramidEst{1});
        subplot(2,3,2); show_flow_as('rgb', f1a.flowPyramidEst{1});
        subplot(2,3,3); show_flow_as('rgb', f1b.flowPyramidEst{1});
        subplot(2,3,4); show_flow_as('rgb', f1.flowPyramidEst{1}-0.5*f1a.flowPyramidEst{1});
        subplot(2,3,5); show_flow_as('rgb', f1.flowPyramidEst{1}-0.5*f1b.flowPyramidEst{1});
        subplot(2,3,6); show_flow_as('rgb', f1a.flowPyramidEst{1}-f1b.flowPyramidEst{1});
    end
end

save(results_path, ...
            'mean_flow_erros',...
            'mean_flow_velocity_erros',...
            'mean_flow_direction_erros',...
            'var_flow_velocity_erros',...
            'var_flow_direction_erros');
%%
figure;
subplot(1,3,1);
hist(abs(mean_flow_erros(:,1)), 100);
title('Mean flow errors - F_0 v F_{1a}');
subplot(1,3,2);
hist(abs(mean_flow_erros(:,2)), 100);
title('Mean flow errors - F_0 v F_{1b}');
subplot(1,3,3);
hist(abs(mean_flow_erros(:,3)), 100);
title('Mean flow errors - F_{1a} v F_{1b}');


figure;
subplot(1,3,1);
hist(mean_flow_velocity_erros(:,1), 100);
title('Mean flow errors (velocity - mean) - F_0 v F_{1a}');
subplot(1,3,2);
hist(mean_flow_velocity_erros(:,2), 100);
title('Mean flow errors (velocity - mean) - F_0 v F_{1b}');
subplot(1,3,3);
hist(mean_flow_velocity_erros(:,3), 100);
title('Mean flow errors (velocity - mean) - F_{1a} v F_{1b}');

figure;
subplot(1,3,1);
hist(mean_flow_direction_erros(:,1), 100);
title('Mean flow errors (direction - mean) - F_0 v F_{1a}');
subplot(1,3,2);
hist(mean_flow_direction_erros(:,2), 100);
title('Mean flow errors (direction - mean) - F_0 v F_{1b}');
subplot(1,3,3);
hist(mean_flow_direction_erros(:,3), 100);
title('Mean flow errors (direction - mean) - F_{1a} v F_{1b}');

figure;
subplot(1,3,1);
hist(var_flow_velocity_erros(:,1), 100);
title('Mean flow errors (velocity - var) - F_0 v F_{1a}');
subplot(1,3,2);
hist(var_flow_velocity_erros(:,2), 100);
title('Mean flow errors (velocity - var) - F_0 v F_{1b}');
subplot(1,3,3);
hist(var_flow_velocity_erros(:,3), 100);
title('Mean flow errors (velocity - var) - F_{1a} v F_{1b}');

figure;
subplot(1,3,1);
hist(var_flow_direction_erros(:,1), 100);
title('Mean flow errors (direction - var) - F_0 v F_{1a}');
subplot(1,3,2);
hist(var_flow_direction_erros(:,2), 100);
title('Mean flow errors (direction - var) - F_0 v F_{1b}');
subplot(1,3,3);
hist(var_flow_direction_erros(:,3), 100);
title('Mean flow errors (direction - var) - F_{1a} v F_{1b}');