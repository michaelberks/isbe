mammo_roi = norm_roi1;
mammo_prob = norm_roi1_prob;
%
[icp_full] = mb_full_image_icp(mammo_roi, 1, 6);
%
[max_icp] = max(icp_full(:,:,:,:), [], 4);
[max_icp] = max(max_icp, [], 3);
%
figure; 
subplot(1,2,1); image(complex2rgb(abs(max_icp) .* exp(2*i*angle(max_icp)))); axis image;
subplot(1,2,2); image(complex2rgb((1-mammo_prob) .* exp(2*i*angle(max_icp)))); axis image;

for level = 1:6
        
    [max_icp] = max(icp_full(:,:,:,level), [], 3);

    subplot(2,3,level);
    image(complex2rgb((1-mammo_prob) .* exp(2*i*angle(max_icp)))); axis image;
    
end
%%
for level = 1:6
        
    [max_icp] = max(icp_full(:,:,:,level), [], 3);

    subplot(2,3,level);
    image(complex2rgb(abs(max_icp) .* exp(2*i*angle(max_icp)))); axis image;
    
end
%%
% [max_icp] = max(icp_full(:,:,:,:), [], 4);
% [max_icp] = max(max_icp, [], 3);