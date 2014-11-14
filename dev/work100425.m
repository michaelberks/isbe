mass_list = dir([mberksroot, 'image_data/masses512x512/*.mat']);

for ii = 1:length(mass_list)
    mass = u_load([mberksroot, 'image_data/masses512x512/' mass_list(ii).name]);
    ff = find(~any(mass), 1, 'first');
    ll = find(~any(mass), 1, 'last');
    
    if ll < 256
%         figure; 
%         subplot(1,2,1); imagesc(mass); axis image; colormap(gray(256));
        mass = padarray(mass(:,ll+5:end), [0 ll+4], 'symmetric', 'pre');
        save([mberksroot, 'image_data/masses512x512/' mass_list(ii).name], 'mass');
%         subplot(1,2,2); imagesc(mass); axis image; colormap(gray(256));
%         size(mass)
    elseif ff > 256
%         figure; 
%         subplot(1,2,1); imagesc(mass); axis image; colormap(gray(256));
        mass = padarray(mass(:,1:ff-5), [0 512-ff+5], 'symmetric', 'post');
        save([mberksroot, 'image_data/masses512x512/' mass_list(ii).name], 'mass');
%         subplot(1,2,2); imagesc(mass); axis image; colormap(gray(256));
%         size(mass)
    end
end
%%
test_ratings = [ones(33,1); 2*ones(6,1); 3*ones(6,1); 4*ones(11,1); 5*ones(2,1);...
                ones(3,1); 2*ones(2,1); 3*ones(2,1); 4*ones(11,1); 5*ones(33,1)];
test_labels = [zeros(58,1); ones(51, 1)];

[roc_pts auc tp_count fp_count auc_se] = calculate_roc_curve(test_ratings,test_labels,0:5);
        