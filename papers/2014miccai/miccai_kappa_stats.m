function [km12 km13 km23] = miccai_kappa_stats(co_occurrence_t, num_boots, sz)

if ~exist('sz', 'var')
    sz = 0;
end

num_images = size(co_occurrence_t,4);

km12_b = zeros(num_boots(1),1);
km13_b = zeros(num_boots(1),1);
km23_b = zeros(num_boots(1),1);

kl12_b = zeros(num_boots(1),1);
kl13_b = zeros(num_boots(1),1);
kl23_b = zeros(num_boots(1),1);

kd12_b = zeros(num_boots(1),1);
kd13_b = zeros(num_boots(1),1);
kd23_b = zeros(num_boots(1),1);
for i_b = 1:num_boots(1)
    if num_boots(1) == 1 && length(num_boots) > 1
        b_idx = 1:num_images;
    else
        b_idx = ceil(num_images*rand(num_images,1));
    end
    co_occurrence_i = co_occurrence_t(:,:,:,b_idx);
    co_occurrence = sum(co_occurrence_i,4);
    co_occurrence(1,1,1) = sz;

    [km12_b(i_b) kl12_b(i_b)] = do_kappa(1,2,3, co_occurrence);
    [km13_b(i_b) kl13_b(i_b)] = do_kappa(1,3,2, co_occurrence);
    [km23_b(i_b) kl23_b(i_b)] = do_kappa(2,3,1, co_occurrence);

    kd12_b(i_b) = do_kappa2(1,2,3, co_occurrence(2:3,2:3,2:3));
    kd13_b(i_b) = do_kappa2(1,3,2, co_occurrence(2:3,2:3,2:3));
    kd23_b(i_b) = do_kappa2(2,3,1, co_occurrence(2:3,2:3,2:3));
    
end

km12 = mean(km12_b);
km13 = mean(km13_b);
km23 = mean(km23_b);

kl12 = mean(kl12_b);
kl13 = mean(kl13_b);
kl23 = mean(kl23_b);

kd12 = mean(kd12_b);
kd13 = mean(kd13_b);
kd23 = mean(kd23_b);

display('--------------------------------------');
display('Agreement on marking a vessel or not')
display(['Agreement between 1 and 2, k_m = ' num2str(km12) ' +\- ' num2str(2*std(km12_b))]);
display(['Agreement between 1 and 3, k_m = ' num2str(km13) ' +\- ' num2str(2*std(km13_b))]);
display(['Agreement between 2 and 3, k_m = ' num2str(km23) ' +\- ' num2str(2*std(km23_b))]);

display(['Difference in agreement between 3 and 1/2, k_m = ' num2str(km12 - (km13 + km23)/2)...
    ' +\- ' num2str(2*std(km12_b - (km13_b + km23_b)/2))]);

display(['Limit of agreement between 1 and 2, k_l = ' num2str(mean(kl12_b)) ' +\- ' num2str(2*std(kl12_b))]);
display(['Limit of agreement between 1 and 3, k_l = ' num2str(mean(kl13_b)) ' +\- ' num2str(2*std(kl13_b))]);
display(['Limit of agreement between 2 and 3, k_l = ' num2str(mean(kl23_b)) ' +\- ' num2str(2*std(kl23_b))]);

display(['Difference in limit agreements between 3 and 1/2, k_m = ' num2str(kl12 - (kl13 + kl23)/2)...
    ' +\- ' num2str(2*std(kl12_b - (kl13_b + kl23_b)/2))]);

display('---');
display('Agreement selecting marked vessels as distal or non-distal')
display(['Agreement between 1 and 2, k_d = ' num2str(kd12) ' +\- ' num2str(2*std(kd12_b))]);
display(['Agreement between 1 and 3, k_d = ' num2str(kd13) ' +\- ' num2str(2*std(kd13_b))]);
display(['Agreement between 2 and 3, k_d = ' num2str(kd23) ' +\- ' num2str(2*std(kd23_b))]);

display(['Difference in agreement between 3 and 1/2, k_d = ' num2str(kd12 - (kd13 + kd23)/2)...
    ' +\- ' num2str(2*std(kd12_b - (kd13_b + kd23_b)/2))]);


function [k k_lim] = do_kappa(ii, jj, kk, co_occurrence)

tc = sum(co_occurrence(:));
oi = squeeze(sum(sum(co_occurrence,jj),kk));
oj = squeeze(sum(sum(co_occurrence,ii),kk));
oij = squeeze(sum(co_occurrence,kk));

mask = logical([1 0 0; 0 1 1; 0 1 1]);
pr_e = (oi(1) / tc) * (oj(1) / tc) + (sum(oi(2:3))/tc * sum(oj(2:3))/tc);
pr_a = sum(oij(mask)) / tc;

k = (pr_a - pr_e) / (1 - pr_e);

b = oij(1,2) + oij(1,3);
c = oij(2,1) + oij(3,1);
d = oij(2,2) + oij(2,3) + oij(3,2) + oij(3,3);

k_lim = 2*d/(2*d + b + c);

function k = do_kappa2(ii, jj, kk, co_occurrence)

tc = sum(co_occurrence(:));
oi = squeeze(sum(sum(co_occurrence,jj),kk));
oj = squeeze(sum(sum(co_occurrence,ii),kk));
oij = squeeze(sum(co_occurrence,kk));

pr_e = (oi(1) / tc) * (oj(1) / tc) + (oi(2) / tc) * (oj(2) / tc);
pr_a = (oij(1,1)+oij(2,2)) / tc;

k = (pr_a - pr_e) / (1 - pr_e);

