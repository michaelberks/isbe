mam = (imread('C:\isbe\density\mammograms\018LCC1824.tif'));
mam = imresize(mam, 0.25, 'bilinear');
step_profile = improfile;
%%
step_profile_mf = medfilt1(step_profile, 25);
step_profile_grad = gradient(step_profile_mf);


for jj = 60:70, 
    figure; plot(step_profile_mf); hold on; 
    for ii = -10:10; 
        plot([2002-ii*jj 2002-ii*jj], [0 2^16], 'r:'); 
    end
end

for jj = 60:70, 
    figure; plot(step_profile_grad); hold on; 
    for ii = -10:10; 
        plot([2002-ii*jj 2002-ii*jj], [-1500 1500], 'r:'); 
    end
end
%%
step_profile_g0 = step_profile_grad;
step_profile_g0(step_profile_grad > 0) = 0;

idx = zeros(length(step_profile_g0), 1);

dd = 20;
for ii = dd+1:length(step_profile_g0)-dd
    idx(ii) = step_profile_g0(ii) == min(step_profile_g0(ii-dd:ii+dd));
end
step_profile_g0(~idx) = 0;

figure; plot(step_profile_g0, 'r+');
%%
for jj = 60:70, 
    figure; plot(step_profile_g0, 'b+'); hold on;
    for ii = -10:10; 
        plot([2002-ii*jj 2002-ii*jj], [-1500 0], 'r:'); 
    end
end
%%
mam = (imread('C:\isbe\density\mammograms\018LCC1824.tif'));
mam = imresize(mam, 0.25, 'bilinear');
figure; imagesc(mam); axis image; colormap(gray(256));
h = imline(gca,[]);
roi = mam(100:300, 200:1000);
figure; imagesc(roi); axis image; colormap(gray(256));
h = impoly(gca,[]);