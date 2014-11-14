%--------------------------------------------------------------------------
% --------------------- Plotting data on images ---------------------------
%--------------------------------------------------------------------------
%%
%Let's look again at the face points we used in the first week, now with an
%actual face

data_dir = 'P:\MatlabTutorial\data\';

%Load in the data
face = imread([data_dir 'face_1.jpg']);
load([data_dir 'face_1.mat'], 'face_pts');

%Display the face image
figure; image(face); axis image; hold on; %What do axis image and hold on do?

%Now, with hold on, we can plot the face pts on the image
plot(face_pts(:,1), face_pts(:,2), 'r-');
%%
%Now try plotting the points on their own
figure; plot(face_pts(:,1), face_pts(:,2), 'r-');

%What do you notice?

%Why do you think Matlab behaves this way?

% -------------------------------------------------------------------------
%% Another example... with mammograms
load([data_dir '024RCC.mat']); %loads in a data structure mammogram
load([data_dir 'an04_024RCC.mat'], 'mass');

%display the image
figure; imagesc(mammogram); axis image; colormap(gray(256)); hold all;

%Mass is a structure (remember them) containing annotations made by a
%radiologist. They were made on an image 4 times the size of this mammogram
%so we need to adjust the coordinates first when we plot. The coordinates
%are also saved with respect to a region around the mass, specified by R1
%and C1 - the first row and column of the RoI respectively

%First plot the mass outline
plot(...
    (mass.mass_outline(:,1) + mass.C1) / 2, ... %adds the column offset to the x-coords and reduces by 2
    (mass.mass_outline(:,2) + mass.R1) / 2, 'r'); %add the row offset to the y-coords

%We also have spicules annotated (a type of structure associated with
%malignant cancers)
for i_sp = 1:length(mass.mass_spicules)
    plot(...
        (mass.mass_spicules(i_sp).outline(:,1) + mass.C1) / 2,...
        (mass.mass_spicules(i_sp).outline(:,2) + mass.R1) / 2);
end
%%
%Finally, we can use axis again, to zoom in on the region just around the
%mass
axis([mass.C1 mass.C2 mass.R1 mass.R2] / 2);

%Add we can use caxis to adjust the contrast, just for this RoI
caxis([min(mass.mass_ROI(:)) max(mass.mass_ROI(:))])

%We can even add text if we want
mass_centre_x = mean((mass.mass_outline(:,1) + mass.C1) / 2);
mass_centre_y = mean((mass.mass_outline(:,2) + mass.R1) / 2);
text(mass_centre_x, mass_centre_y, 'A malignant cancer');






