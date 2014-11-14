for sigma = [1 2 4]
    display(['Sigma =' num2str(sigma)]);

    [g dg] = gaussian_filters_1d(sigma);
    x = 1:100;

    for k = [pi/4 0.5 0.2 0.1 0.01]
        theta = x*k;
        ctheta = complex(cos(2*theta), sin(2*theta));

        k_est = real( -.5i*conj(ctheta).*conv2(ctheta, dg/sigma, 'same') );
        display(['k = ' num2str(k) ', k_est = ' num2str(k_est(50))])
    end
end
%%
x = repmat(-64:63, 128, 1);
y = x';

theta_xy = -atan2(y, x);
figure; imgray(theta_xy); colormap(hsv(256));

I_xy = cos(sqrt(x.^2 + y.^2));
figure; imgray(I_xy);

k_est = complex_gaussian_curvature(theta_xy, 1);
figure; imgray(k_est);
figure; imgray(cos(1./k_est));