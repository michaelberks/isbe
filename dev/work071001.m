%%
%
% Make ratio errors of the scale model loo erros
%

load C:\isbe\dev\scale\errors

ratio_indie_s = ers_indie2(:,1)./ ers_indie2(:,2);
ratio_shape_s = ers_shape2(:,1)./ ers_shape2(:,2);
ratio_tex_s = ers_tex2(:,1)./ ers_tex2(:,2);
ratio_scale_s = ers_scale2(:,1)./ ers_scale2(:,2);
ratio_combined_s = ers_combined2(:,1)./ ers_combined2(:,2);

ratio_indie_l = ers_indie3(:,1)./ ers_indie3(:,2);
ratio_shape_l = ers_shape3(:,1)./ ers_shape3(:,2);
ratio_tex_l = ers_tex3(:,1)./ ers_tex3(:,2);
ratio_scale_l = ers_scale3(:,1)./ ers_scale3(:,2);
ratio_combined_l = ers_combined3(:,1)./ ers_combined3(:,2);

ratio_indie_r50 = ers_indie4(:,1)./ ers_indie4(:,2);
ratio_shape_r50 = ers_shape4(:,1)./ ers_shape4(:,2);
ratio_tex_r50 = ers_tex4(:,1)./ ers_tex4(:,2);
ratio_scale_r50 = ers_scale4(:,1)./ ers_scale4(:,2);
ratio_combined_r50 = ers_combined4(:,1)./ ers_combined4(:,2);

ratio_indie_r51 = ers_indie5(:,1)./ ers_indie5(:,2);
ratio_shape_r51 = ers_shape5(:,1)./ ers_shape5(:,2);
ratio_tex_r51 = ers_tex5(:,1)./ ers_tex5(:,2);
ratio_scale_r51 = ers_scale5(:,1)./ ers_scale5(:,2);
ratio_combined_r51 = ers_combined5(:,1)./ ers_combined5(:,2);

%save the error matrices
save C:\isbe\dev\scale\ratios ratio*
%%
%
% Plot the results in boxplots
%
ratio_indie = nan; ratio_indie = ratio_indie(ones(51, 4));
ratio_shape = ratio_indie;
ratio_tex = ratio_indie;
ratio_scale = ratio_indie;
ratio_combined = ratio_indie;

% individual errors boxplot
figure;
ratio_indie(1:50, [1 3]) = [ratio_indie_s, ratio_indie_r50];
ratio_indie(:, [2 4]) = [ratio_indie_l, ratio_indie_r51];
boxplot(ratio_indie);
title('individual errors boxplot');

% shape errors boxplot
figure;
ratio_shape(1:50, [1 3]) = [ratio_shape_s, ratio_shape_r50];
ratio_shape(:, [2 4]) = [ratio_shape_l, ratio_shape_r51];
boxplot(ratio_shape);
title('shape errors boxplot');

% tex errors boxplot
figure;
ratio_tex(1:50, [1 3]) = [ratio_tex_s, ratio_tex_r50];
ratio_tex(:, [2 4]) = [ratio_tex_l, ratio_tex_r51];
boxplot(ratio_tex);
title('tex errors boxplot');

% scale errors boxplot
figure;
ratio_scale(1:50, [1 3]) = [ratio_scale_s, ratio_scale_r50];
ratio_scale(:, [2 4]) = [ratio_scale_l, ratio_scale_r51];
boxplot(ratio_scale);
title('scale errors boxplot');

% combined errors boxplot
figure;
ratio_combined(1:50, [1 3]) = [ratio_combined_s, ratio_combined_r50];
ratio_combined(:, [2 4]) = [ratio_combined_l, ratio_combined_r51];
boxplot(ratio_combined);
title('combined errors boxplot');

%%
% Look at the mean ratio errors
[mean(ratio_indie_s), mean(ratio_indie_l), mean(ratio_indie_r50), mean(ratio_indie_r51)]
[mean(ratio_shape_s), mean(ratio_shape_l), mean(ratio_shape_r50), mean(ratio_shape_r51)]
[mean(ratio_tex_s), mean(ratio_tex_l), mean(ratio_tex_r50), mean(ratio_tex_r51)]
[mean(ratio_scale_s), mean(ratio_scale_l), mean(ratio_scale_r50), mean(ratio_scale_r51)]
[mean(ratio_combined_s), mean(ratio_combined_l), mean(ratio_combined_r50), mean(ratio_combined_r51)]