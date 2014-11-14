function [] = mb_change_mass_region_size(varargin)
%MB_CHANGE_MASS_REGION_SIZE *Insert a one line summary here*
%   [] = mb_change_mass_region_size(varargin)
%
% MB_CHANGE_MASS_REGION_SIZE uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%
% Optional Arguments:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 28-May-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin, '0',...
    {... %mandatory arguments
    'MassList',...
    'SavePath'...
    },... %optional argument and defaults
    'AnnoPath', 'C:\isbe\dev\annotations\',...
    'Size', [1024 1024],... 
    'AnnoPrefix', 'an',...
    'Pad', 1,...
    'Plot', 0 ...
    );

%make sure directory paths are all filesep ended
if ~strcmp(args.AnnoPath(end), filesep);
    args.AnnoPath = [args.AnnoPath, filesep];
end


for ii = 1:length(args.MassList)
    %Get the annotated mass name from the mass file list name (i.e. take
    %off the 'm' from the start of the mass name and add the annotation
    %prefix)
    mass_name = [args.AnnoPath, args.AnnoPrefix, args.MassList(ii).name(2:end)];
    
    %load in the annotated mass
    mass_anno = u_load(mass_name);
    
    %load in the mamogram from which the masses was taken (assuming here
    %that this is an image and so use imread)
    mammo = double(imread(mass_anno.name));
    [r c] = size(mammo);
    
    %Get the geometric centre of the annotated mass, relative to the
    %mammogram
    gc = mean(mass_anno.mass_outline) + [mass_anno.C1 mass_anno.R1] - 1;
    
    %Get start and end of new bounding box
    r1 = max(1, round(gc(2) - args.Size(1)/2));
    c1 = max(1, round(gc(1) - args.Size(2)/2));
    r2 = min(r, round(gc(2) + args.Size(1)/2)-1);
    c2 = min(c, round(gc(1) + args.Size(2)/2)-1);
    
    %Now start building the new mass structure:
    %First get the new mass ROI
    mass = [];
    if args.Pad
        mass.mass_ROI = zeros(args.Size);
        mass.mass_ROI(1:r2-r1+1, 1:c2-c1+1) = mammo(r1:r2, c1:c2);
    else
        mass.mass_ROI = mammo(r1:r2, c1:c2);
    end
    
    %Then import the mass background into the mammogram, and get the mass
    %background
    mammo(mass_anno.R1:mass_anno.R2, mass_anno.C1:mass_anno.C2) = ...
        mammo(mass_anno.R1:mass_anno.R2, mass_anno.C1:mass_anno.C2) - ...
        mass_anno.mass_sub_it;
    
    mass.background_ROI = mammo(r1:r2, c1:c2);
    
    %Can clear the mammo from memmory now
    clear mammo;
    
    %Now we need to update the mass_outline to the new coordinate
    %system of the region
    mass.mass_outline = mass_anno.mass_outline;
    mass.mass_outline(:,1) = mass.mass_outline(:,1) + mass_anno.C1 - c1;
    mass.mass_outline(:,2) = mass.mass_outline(:,2) + mass_anno.R1 - r1;
    
    %Do the same for any mass spicules
    mass.mass_spicules = mass_anno.mass_spicules;
    for jj = 1:length(mass.mass_spicules)
        mass.mass_spicules(jj).outline(:,1) = ...
            mass.mass_spicules(jj).outline(:,1) + mass_anno.C1 - c1;
        mass.mass_spicules(jj).outline(:,2) = ...
            mass.mass_spicules(jj).outline(:,2) + mass_anno.R1 - r1;
    end
    
    %Finally, flip the new mass if it is from the right breast
    if ~isempty(strfind(args.MassList(ii).name, 'R'))
        mass.mass_outline(:,1) = args.Size(2) + 1 - mass.mass_outline(:,1);
        mass.mass_ROI = fliplr(mass.mass_ROI);
        mass.background_ROI = fliplr(mass.background_ROI);
        
        for jj = 1:length(mass.mass_spicules)
            mass.mass_spicules(jj).outline(:,1) = ...
                args.Size(2) + 1 - mass.mass_spicules(jj).outline(:,1);
        end
        %display(['Flipping ', args.MassList(ii).name]);
    end
    
    if args.Plot
        figure;
        subplot(1,2,1); imagesc(mass.mass_ROI); axis image; colormap(gray(256));
        hold on; plot(mass.mass_outline(:,1), mass.mass_outline(:,2));
        subplot(1,2,2); imagesc(mass.background_ROI); axis image; colormap(gray(256));
        hold on; plot(mass.mass_outline(:,1), mass.mass_outline(:,2));
    end
    
    if ~isempty(args.SavePath)
        %make sure directory paths are all filesep ended
        if ~strcmp(args.SavePath(end), filesep);
            args.SavePath = [args.SavePath, filesep];
        end
        %Save the mass to it's new home
        save([args.SavePath, args.MassList(ii).name], 'mass');
    end
end
    
    
    
    
    

