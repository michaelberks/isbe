%Make test data
mass_list = dir('C:\isbe\dev\annotations\*.mat');
normal_list = dir('C:\isbe\dev\background\images\normal1024\*.bmp');
r_mass = randperm(length(mass_list));
r_normal = randperm(length(normal_list));

%
mkdir C:\isbe\thesis\figures\5\experiment

mm = 1; nn = 1;
while nn <= 50
    bg = imread(['C:\isbe\dev\background\images\normal1024\', normal_list(r_normal(nn)).name]);
    mass = u_load(['C:\isbe\dev\annotations\', mass_list(r_mass(mm)).name]);
    [r c] = size(mass.mass_ROI);
    if r <= 1024 && c < 1024
        bg = double(bg(1:r, 1:c));
        test_data.mass_ROI = bg + mass.mass_sub_it;
        test_data.mass_outline = mass.mass_outline;
        mass_bw = roipoly(bg, mass.mass_outline(:,1), mass.mass_outline(:,2));
        test_data.real_bg = bg(mass_bw);
        save(['C:\isbe\thesis\figures\5\experiment\test_mass', zerostr(nn,2)], 'test_data');
        clear test_data
        nn = nn+1;
    end
    
    mm = mm+1;
end


%%
test_list = dir('C:\isbe\thesis\figures\5\experiment\*.mat');

[n1 sigma] = meshgrid(0:10:100, 5:5:30);
errors_original = zeros(size(n1));
errors_iterative = zeros(size(n1));
%%
for ii = 37:numel(n1)
    [errors] = subtract_mass_test(test_list,'C:\isbe\thesis\figures\5\experiment\', n1(ii), 20, 10, sigma(ii), 1, 'biharmTPS');

    errors_original(ii) = mean(errors(:,1));
    errors_iterative(ii) = mean(errors(:,2));
    save C:\isbe\thesis\figures\5\experiment\errors errors_*
end
%%
for ii = 1:11
    [errors] = subtract_mass_test(test_list,'C:\isbe\thesis\figures\5\experiment\', n1(1,ii), 20, 10, 35, 1, 'biharmTPS');

    errors_original(7,ii) = mean(errors(:,1));
    errors_iterative(7,ii) = mean(errors(:,2));
    save C:\isbe\thesis\figures\5\experiment\errors errors_*
end
%%
n2 = 10:10:50;
for ii = 1:length(n2)
    [errors] = subtract_mass_test(test_list,'C:\isbe\thesis\figures\5\experiment\', 60, n2(ii), 10, 25, 1, 'biharmTPS');
    errors_n2(ii) = mean(errors(:,1));
end

%%
d = 10:10:50;
for ii = 1:length(d)
    [errors] = subtract_mass_test(test_list,'C:\isbe\thesis\figures\5\experiment\', 60, 20, d(ii), 25, 1, 'biharmTPS');
    errors_d(ii) = mean(errors(:,1));
end
%%
%original best
[min_original i_o] = min(errors_original(:));
display(['Best original parameters are n1 = ', num2str(n1(i_o)), ' sigma = ', num2str(sigma(i_o))]);

%best iterative is:
[min_iterative i_i] = min(errors_iterative(:));
display(['Best iterative parameters are n1 = ', num2str(n1(i_i)), ' sigma = ', num2str(sigma(i_i))]);