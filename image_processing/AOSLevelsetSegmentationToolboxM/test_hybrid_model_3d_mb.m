%%
subject_dir = 'Q:/data/MB/travastin/A011_CJ/';
tumour_path = 'FA_1/tumour/Tofts_plus_vp_Auto_AIF_results/';

visit = 1;
visit_dir = [subject_dir 'visit' num2str(visit) '/'];
dynamic_dir = [visit_dir 'dynamic/'];


[dyn_vol, dyn_header] = load_img_volume([dynamic_dir 'dyn_50.hdr']);
tumour_mask = load_img_volume([visit_dir tumour_path 'ROI.raw.hdr']) > 0;

x_proj = any(any(tumour_mask,3),1);
y_proj = any(any(tumour_mask,3),2);
z_proj = squeeze(any(any(tumour_mask,2),1));

x_idx = (find(x_proj, 1, 'first')-2):(find(x_proj, 1, 'last')+2);
y_idx = (find(y_proj, 1, 'first')-2):(find(y_proj, 1, 'last')+2);
z_idx = (find(z_proj, 1, 'first')-1):(find(z_proj, 1, 'last')+1);

%V = dyn_vol(y_idx, x_idx, z_idx); 
%phi = double(tumour_mask(y_idx, x_idx, z_idx));

V = max(dyn_vol(:)) - dyn_vol; 
phi = zeros(size(V));
phi(y_idx, x_idx, z_idx) = 1;
%%
propagation_weight = 1e-3; 
GAC_weight = .02; 
% g = ones(size(V)); % linear diffusion 
g = ac_gradient_map(V,1,0,5*ones(1,3), 2*eye(3)); 
delta_t = 4; 
mu = 1200; 
%%
for i = 1:10    
    phi = ac_hybrid_model(V-mu, phi-.5, propagation_weight, GAC_weight, g, ...
        delta_t, 1); 
    figure;
    iso = isosurface(phi,0);
    h = patch(iso,'edgecolor','r','facecolor','w');  axis equal;  view(3); 
    set(gcf,'name', sprintf('#iters = %d',i));
    drawnow; 
end

%%
figure;
slice = 10:17;
for i = 1:8
    subplot(2,4,i); imshow(V(:,:,slice(i)),[]); hold on; 
    c = contours(phi(:,:,slice(i)),[0,0]);
    zy_plot_contours(c,'linewidth',2);
end