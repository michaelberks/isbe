theta_pi = mod(theta, pi);
for my_ang = 0:12
    my_oris = pi*(my_ang/12)*ones(181);
    theta_diff = abs(theta_pi-my_oris);
    theta_diff(rij < r_min) = inf;
    figure;
    subplot(2,2,1); imagesc(theta_diff); axis image; colorbar;
    subplot(2,2,2); imagesc(abs(theta_diff-pi)); axis image; colorbar;
    subplot(2,2,3); imagesc(theta_diff < phi); axis image;
    subplot(2,2,4); imagesc(abs(theta_diff-pi) < phi); axis image;
end
    