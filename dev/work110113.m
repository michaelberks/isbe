[orientation_errors line_contrasts] = compute_image_orientation_errors(...
    test_dir, [prob_dir param_dir], 'centre_line');
%%
con_min = floor(min(line_contrasts));
con_max = floor(max(line_contrasts));

mean_errors = zeros(con_max-con_min+1,1);
for ii = con_min:con_max
    idx = (line_contrasts >= ii) & (line_contrasts < ii+1);
    mean_errors(ii) = mean(abs(orientation_errors(idx)));
end
%%
[sorted_contrasts idx] = sort(line_contrasts);
sorted_errors = abs(orientation_errors(idx));
%%
smoothed_errors = imfilter(sorted_errors, ones(1e4,1)/1e4, 'same', 'replicate');

 plot(sorted_contrasts, smoothed_errors, 'r.'); hold on;
%%
test_dir = 'C:\isbe\asymmetry_project\data\synthetic_lines\real512\';
prob_dir = 'C:\isbe\asymmetry_project\data\synthetic_lines\real512\results\';
%
figure; hold on;
%
[orientation_errors line_contrasts] = compute_image_orientation_errors(...
    test_dir, [prob_dir '233902'], 'centre_line');
save('C:\isbe\asymmetry_project\data\misc\ori_errors_233902.mat', 'orientation_errors');
save('C:\isbe\asymmetry_project\data\misc\line_contrasts.mat', 'line_contrasts');

[sorted_contrasts idx] = sort(line_contrasts);
sorted_errors = abs(orientation_errors(idx));
smoothed_errors = imfilter(sorted_errors, ones(1e4,1)/1e4, 'same', 'symmetric');
plot(sorted_contrasts, smoothed_errors, 'r');
%%
[orientation_errors] = compute_image_orientation_errors(...
    test_dir, [prob_dir '233908'], 'centre_line');
[sorted_contrasts idx] = sort(line_contrasts);
sorted_errors = abs(orientation_errors(idx));
smoothed_errors = imfilter(sorted_errors, ones(1e4,1)/1e4, 'same', 'symmetric');
plot(sorted_contrasts, smoothed_errors, 'g');
%
[orientation_errors] = compute_image_orientation_errors(...
    test_dir, [prob_dir 'g2d_scales_16\orientations'], 'centre_line');
[sorted_contrasts idx] = sort(line_contrasts);
sorted_errors = abs(orientation_errors(idx));
smoothed_errors = imfilter(sorted_errors, ones(1e4,1)/1e4, 'same', 'symmetric');
plot(sorted_contrasts, smoothed_errors, 'b');