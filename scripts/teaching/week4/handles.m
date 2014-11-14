%--------------------------------------------------------------------------
% --------------------- Handles ---------------------------
%--------------------------------------------------------------------------
%%
%Let's repeat the mammogram figure, but we'll set the display properties of
%each figure

%% Mammogram figure
data_dir = 'P:\MatlabTutorial\data\';
load([data_dir '024RCC.mat']); %loads in a data structure mammogram
load([data_dir 'an04_024RCC.mat'], 'mass');
%%
%display the image
f1 = figure; %Creates a new figure object, and save its handle into f1
a1 = axes; %Create a new axes object in f1, saves the handle in f2
i1 = image(mammogram); %Creates a new image object in a1

%before continuing, lets look at the properties we can set for each object
display('*********** Figure properties *******************');
set(f1)
display('*********** Axes properties *******************');
set(a1)
display('*********** Image properties *******************');
set(i1)
%%
%Set the figure's colormap to gray
set(f1, 'colormap', gray(256));

%Now lets set a bunch of stuff in the axes
set(a1, ...
    'YDir', 'reverse', ... %This inverts the y-axis, like axis image (actually, this will already be set, as we used 'image')
    'DataAspectRatio', [1 1 1], ... %Make the axis equal scale
    'NextPlot', 'add'... %Switches hold on
    );
setappdata(a1,'PlotHoldStyle',true); %This causes different colours to be used for each line (like hold all)

%We can now plot the mass outline and spicules
p1 = plot(...
    (mass.mass_outline(:,1) + mass.C1) / 2, ... %adds the column offset to the x-coords and reduces by 4
    (mass.mass_outline(:,2) + mass.R1) / 2, 'r'); %add the row offset to the y-coords
for i_sp = 1:length(mass.mass_spicules)
    plot(...
        (mass.mass_spicules(i_sp).outline(:,1) + mass.C1) / 2,...
        (mass.mass_spicules(i_sp).outline(:,2) + mass.R1) / 2);
end
%%
%Finally, to zoom in round the mass, and rescale the colour map
set(a1,...
    'Xlim', [mass.C1 mass.C2] / 2,...
    'Ylim', [mass.R1 mass.R2] / 2,...
    'Clim', [min(mass.mass_ROI(:)) max(mass.mass_ROI(:))]);
%%
%Who belongs to who? Use children/parent properties
display([a1 get(f1, 'children')]); %The axes belongs to the figure
display([f1 get(a1, 'parent')]);
display([a1 get(i1, 'parent')]); %The image belongs to the axes
display([a1 get(p1, 'parent')]); %As does the plot of the mass outline
display(get(a1, 'children')); %The axes has lots of other children too - all the spicules!

%% Ok, so far so good, but if all we're doing is replicating things we can do with shortcut functions, what's the point?

%What if we want to change, where the axes are positioned in the image
set(a1, 'position', [0 0 1 1]);

%Or change the order of colours that get used for each new line
set(a1, 'colororder', summer(length(mass.mass_spicules)));
for i_sp = 1:length(mass.mass_spicules)
    plot(...
        (mass.mass_spicules(i_sp).outline(:,1) + mass.C1) / 2,...
        (mass.mass_spicules(i_sp).outline(:,2) + mass.R1) / 2);
end

%Or doing lots of things! 






