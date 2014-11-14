%This script attempts to rotate mass by an arbitrary angle

%Step 1 is extending the region so we don't get a 'zero region' when we
%rotate - this would cause nasty edge effects when we apply the DT-CWT.

load C:\isbe\dev\files\u_files.mat
target_list = dir('C:\isbe\dev\background\images\normal1024\*.bmp');

%load in a test mass and background region - note we have chosen these
%because the dominant orientation of the mass and background textures clash
% see * below for an example in showing how simply reflecting the mass
% vastly improves the appearance of the mass structure transferred into the
% background....
mass = u_load(['C:\isbe\dev\masses\', u_files1(2).name]);
bg = double(imread(['C:\isbe\dev\background\images\normal1024\', target_list(4).name]));

mass_rotate0 = double(mass.mass_ROI) - mass.subtract_ROI;

%To test we're correctly aligning the mass_outline after rotation lets work
%with a mask of the mass
mass_bw = roipoly(mass_rotate0, mass.mass_outline(:,1), mass.mass_outline(:,2));

%Extend the mask, rotate, then take take the central region, the resulting
%region should have the rotated mass centered as in the original region
[r c] = size(mass_rotate0);
pad_vec = ceil(0.5*(sqrt(r^2 + c^2) - [r c]));
%%
degrees = 20;
mass_bw_pad = padarray(mass_bw, pad_vec, 'symmetric');
mass_bw45_pad = imrotate(mass_bw_pad, degrees, 'nearest', 'crop');
mass_bw45 = mass_bw45_pad(pad_vec(1)+(1:r), pad_vec(2)+(1:c));

figure; imagesc(mass_bw); axis image; colormap(gray(256));
figure; imagesc(mass_bw45_pad); axis image; colormap(gray(256));
figure; imagesc(mass_bw45); axis image; colormap(gray(256));

%The mass_outline is defined relative to the top-left corner of the image.
%A rotation through theta using imrotate, w.r.t to the mass_outline
%co-ordinate system, is thus a rotation of -theta about (c/2, r/2). So the
%rotated outline is obtained by (x,y)' = R(-theta)((x,y) - (c/2,r/2)) + (c2,r/2)
theta = pi*degrees / 180;
r_theta = [cos(theta) sin(theta); -sin(theta) cos(theta)];
rot_outline = (r_theta * ...
    (mass.mass_outline - repmat([c r] / 2, size(mass.mass_outline, 1), 1))')' +...
    repmat([c r] / 2, size(mass.mass_outline, 1),1);

hold on;
plot(rot_outline(:,1), rot_outline(:,2), 'r');

%%
%Now try transferring the rotated region - actually let's add this into to
%the transfer function so we don't have to rotate the mass_ROI and the
%subtract_ROI fields. We can also get away with simply rotated the BW mask
%rather than rotating the mass_outline and creating a new mask...

args.MassDilate = 50;
args.Mass = mass;
args.TargetImage = bg;
args.Location = ([1024 1024] - fliplr(size(args.Mass.mass_ROI))) / 2;
[modified_region0] = mb_transfer_mass_dual_tree(args);
args.Rotate = pi/2;
[modified_region45] = mb_transfer_mass_dual_tree(args);

figure; imagesc(modified_region0); axis image; colormap(gray(256));
figure; imagesc(modified_region45); axis image; colormap(gray(256));
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%* Note the clash in orientation between background and mass texture and
%the improved appearance of transferring the mass into the flipped
%background
args.TargetImage = fliplr(bg);
[modified_region_flip] = mb_transfer_mass_dual_tree(args);
figure; imagesc(fliplr(modified_region_flip)); axis image; colormap(gray(256));
