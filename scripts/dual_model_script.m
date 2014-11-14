n_angles = 4;
sigma_x = 2;
sigma_y = sigma_x;
tau = sigma_x*2.35;
halfwidth = 5 * max(sigma_x, sigma_y);
remove_bias = 0;


gabor_f = real(gabor_filters(1, sigma_x, sigma_y, tau, halfwidth, remove_bias));
gauss_f = gaussian_filters_1d(sigma_x) / sigma_x;

% figure;  
% subplot(1,2,1); imgray(gabor_f);
% subplot(1,2,2); imgray(gauss_f);

fg_heights = linspace(1,10,100);
bg_heights = linspace(1,10,100);
k = 0.8;
tau = 2;
gabor_responses = zeros(100);
gauss_responses = zeros(100);
for i_bg = 1:100
    
 
    for i_fg = 1:100
        [line] = ...
            create_ellipse_bar(1, fg_heights(i_fg), 90, 21, 21, 11, 11);
        test_im = line + bg_heights(i_bg);
        gabor_response_i = conv2(test_im, gabor_f, 'valid');
        gauss_response_i = conv2(gauss_f', gauss_f, test_im, 'valid');
        
        gabor_responses(i_bg,i_fg) = gabor_response_i;
        gauss_responses(i_bg,i_fg) = gauss_response_i;
        
    
    end
end
%
figure; 
imagesc(gabor_responses);
title('Gabor responses');
ylabel('Background height');
xlabel('Line height');

figure; 
imagesc(gauss_responses);
title('Gaussian responses');
ylabel('Background height');
xlabel('Line height');  
            
for alpha = linspace(0,1,10)
    Rp = gabor_responses + alpha*gauss_responses;
    Rn = gabor_responses - alpha*gauss_responses;
    sigmoid = 1 ./ (1 + exp(-2*k*Rn));
    dual_model = Rp .* sigmoid;

    figure('Name', ['Models for \alpha = ' num2str(alpha)]);
    subplot(1,3,1);
    imagesc(Rp);
    title('R_p');
    ylabel('Background height');
    xlabel('Line height');
    subplot(1,3,2);
    imagesc(sigmoid);
    title('R_n sigmoid');
    ylabel('Background height');
    xlabel('Line height');
    subplot(1,3,3);
    imagesc(dual_model);
    title('Dual model');
    ylabel('Background height');
    xlabel('Line height');
end
    
% Cp = gabor_responses + tau;
% Cn = gabor_responses - tau;
% constant_model = Cp ./ (1 + exp(-2*k*Cn));
% 
% figure;
% imagesc(Rp);
% title('Constant model responses');
% xlabel('Line height');
% ylabel('Response');

%%
k = 1;
[a b] = meshgrid(linspace(-0,20, 100), linspace(1,10,100));
dm = (a + b) ./ (1 + exp(-2*k*(a-b)));
figure; surf(a, b, dm);
xlabel('Gabor');
ylabel('\alpha*Gaussian');
figure; plot(dm(1:10:end,:)');

figure; surf(a, b, a+b);
figure; surf(a, b, 1./(1 + exp(-2*k*(a-b))));


    