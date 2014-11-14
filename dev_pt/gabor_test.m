close all;

%%
% Gabor filter responses for synthetic lines at different scales - try % swapping over the normalising constant in gabor_filters for f_type = 1:2
for f_type = 1:2
    if f_type == 1
        feat = 'line';
    else
        feat = 'edge';
    end
    
    f1 = figure;
    a1 = subplot(2,2,1); hold all; 
    title(['Abs Gabor response to ' feat]); xlabel('\sigma of filter');
    a2 = subplot(2,2,2); hold all; 
    title(['G" response to ' feat]); xlabel('\sigma of filter');
    a3 = subplot(2,2,3); hold all; 
    title(['Real Gabor response to ' feat]); xlabel('\sigma of filter');
    a4 = subplot(2,2,4); hold all; 
    title(['Imag Gabor response to ' feat]); xlabel('\sigma of filter');
    
    for line_width = 1:8
        if f_type == 1
            f_image = create_sin_bar(line_width, 1, 0, 32, 32,0, 16, 16);
        else
            f_image = create_sin_step(line_width, 1, 0, 32, 32, 16, 16);
        end

%         f_image = -f_image;
        f_image = f_image + 100;
        
        scale_responses_gabor = zeros(8,1);
        scale_responses_gauss = zeros(8,1);
        for sigma = 1:8
            [gabor_responses] = compute_gabor_responses(f_image, sigma, 2);
            [gauss_responses] = gaussian_2nd_derivative_line(f_image, sigma);
            scale_responses_gabor(sigma) = max(gabor_responses(16,16,:,:), [], 4);
            scale_responses_gauss(sigma) = abs(gauss_responses(16,16));
            
%             if line_width == sigma
%                 figure;
%                 %Take maximum response across orientations
%                 [max_responses] = max(gabor_responses(:,:,:,:), [], 4);
%                 %Look at abs, real and imag parts of responses
%                 subplot(1,3,1); imgray(abs(max_responses)); colorbar('southoutside'); 
%                 title(['Absolute response, \sigma = ' num2str(sigma)]);
%                 subplot(1,3,2); imgray(real(max_responses)); colorbar('southoutside');
%                 title(['Real response, \sigma = ' num2str(sigma)]);
%                 subplot(1,3,3); imgray(imag(max_responses)); colorbar('southoutside');
%                 title(['Imaginary response, \sigma = ' num2str(sigma)]);
%             end
        end
        figure(f1);
        plot(a1, 1:8, abs(scale_responses_gabor), '-x');
        plot(a2, 1:8, scale_responses_gauss, '-x');
        plot(a3, 1:8, real(scale_responses_gabor), '-x');
        plot(a4, 1:8, imag(scale_responses_gabor), '-x');
            
    end
    legend(a1, strcat('Feature width =  ', num2str((1:8)')));
    legend(a2, strcat('Feature width =  ', num2str((1:8)')));
    legend(a2, strcat('Feature width =  ', num2str((1:8)')));
    legend(a3, strcat('Feature width =  ', num2str((1:8)'))); 
end 

% return

%% 
%Gabor responses for the retinograms 

%Load retinogram 
ret = u_load([asymmetryroot 'data\retinograms\DRIVE\training\images\21_training.mat']);

% return

newim = zeros(size(ret));
ret = ret(:,:,2);
clear a;
for sigma = [4]
    %Compute Gabor responses
    [gabor_responses] = compute_gabor_responses(ret, sigma, 6);
    
    %Take maximum response across orientations
    [max_responses] = max(gabor_responses(:,:,:,:), [], 4);
    %Look at abs, real and imag parts of responses
    figure;
    a(1) = subplot(2,2,1); imgray(ret);
    title(['Original image']);
    a(2) = subplot(2,2,2); imgray(abs(max_responses)); 
    title(['Absolute response, \sigma = ' num2str(sigma)]);
    a(3) = subplot(2,2,3); imgray(real(max_responses));
    title(['Real response, \sigma = ' num2str(sigma)]);
    a(4) = subplot(2,2,4); imgray(imag(max_responses));
    title(['Imaginary response, \sigma = ' num2str(sigma)]);
    linkaxes(a);
    
    newim(:,:,1) = real(max_responses);
    newim(:,:,2) = imag(max_responses);
    newim(:,:,3) = 0;
    
    newim = newim / max(abs(newim(:)));
    newim = 0.5*newim + 0.5;
    newim(:,:,3) = 0;
    
    figure(); image(newim);
    
    %You might need to zoom to see the responses at vessel centres - the
    %axes are linked
end
