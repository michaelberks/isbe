% SR_demo.m
%
% http://www.robots.ox.ac.uk/~elle/SRcode/
%
% (Just run this -- should be self-explanatory!)
%

%%
% 1: Get a (synthetic) super-res problem:
[o,gtruth] = synthdata_demo(20); % Get struct "o" containing SR data, and "gtruth" ground truth image.
[biv,bih] = size(gtruth);
K = numel(o);
%%
% 2: Look at the data:
figure;
for i = 1:min(K,6)
    subplot(2,3,i); imgray(o(i).im);
    title(['low-res ' num2str(i)]);
end
%%
% 3: Get matrices and set up options vector/params:
[W, Y, La, Lb, M] = makeW(biv,bih,o);
opts = zeros(1,18); % This is the "options" vector for the Netlab "scg" routine.
opts(1) = 1; % verbose
opts(2:3) = 1e-3; % convergeance criteria
opts(14) = 100; % number of iterations before automatic termination.
alp = 0.08;
nu = 0.04;
%%

% 4: Find Average Image, Maximum Likelihood and Huber super-res image estimates:
[avim,msk,M] = getAvim(biv,bih,o);
im_ml = superres_ml(W,Y,La,Lb,avim);%,opts
im_huber = superres_huber(W,Y,La,Lb,avim,alp,nu,opts);

% 5: Look at the various outputs, comparing them to the ground truth image:
figure;
subplot(2,2,1); imgray(avim((gap+1):end-gap,(gap+1):end-gap)+0.5);
title('Average Image');
subplot(2,2,2); imgray(im_ml((gap+1):end-gap,(gap+1):end-gap)+0.5);
title('ML Image');
subplot(2,2,3); imgray(im_huber((gap+1):end-gap,(gap+1):end-gap)+0.5);
title('Huber Image');
subplot(2,2,4); imgray(gtruth((gap+1):end-gap,(gap+1):end-gap)+0.5);
title('Ground Truth Image');


% 6: Be happy :-)
