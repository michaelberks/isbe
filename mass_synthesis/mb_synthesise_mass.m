function [mass rotation_angle] = ...
    mb_synthesise_mass(target_region, target_centre, mass_model, mass_dir, mass_list, mass_idx, rotation_angle)
%MB_SYNTHESISE_MASS *Insert a one line summary here*
%   [mass] = mb_synthesise_mass(region)
%
% Inputs:
%      region- *Insert description of input variable here*
%
%
% Outputs:
%      mass- *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 29-May-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

%This is the final function that given a region and a mass_model

if nargin < 6
    mass_idx = [];
end
if nargin < 7
    rotation_angle = [];
end

%Check size of model matches size of mass_list
if length(mass_list) ~= size(mass_model.B_scale,2)
    error('The number of masses in the list must match the number of masses in the model');
end

%Randomly select a real mass to act as a template to modify the texture in
%the target region
if isempty(mass_idx)
    mass_idx = ceil(rand*length(mass_list));
    display(['Template mass_selected = ', num2str(mass_idx)]);
end

%make sure mass directory ends in filesep
if ~strcmp(mass_dir(end), filesep);
    mass_dir = [mass_dir, filesep];
end
%Load in the template mass
mass_template = u_load([mass_dir, mass_list(mass_idx).name]);

%If we haven't been given a rotation angle we need to work this out...
if isempty(rotation_angle)
    %Compute the orientation map of the template mass
    [template_orientation_map] = ...
        mb_dt_orientation_map('Image', mass_template.background_ROI);
    
    n_bins = 120;
    template_centre = mean(mass_template.mass_outline);
    max_radius = max(abs([...
        mass_template.mass_outline(:,1) - template_centre(1);...
        mass_template.mass_outline(:,2) - template_centre(2)])) + 20;
    
    [x y] = meshgrid(1:size(mass_template.mass_ROI,2),(1:size(mass_template.mass_ROI,1))'); 
    %Compute the orientation hist of the mass region 
    mass_mask = (x-target_centre(1)).^2 + (y-target_centre(2)).^2 < max_radius^2;

    [tout rout template_orientation_hist] = ...
        weighted_complex_rose(template_orientation_map(mass_mask), n_bins);
    clear dummy;
    
%     figure;
%     subplot(1,2,1); polar(tout, rout); title('Template orientation')
    
    %Compute the orientation map of the target region
    [target_orientation_map] = ...
        mb_dt_orientation_map('Image', target_region); 
    
    [x y] = meshgrid(1:size(target_region,2),(1:size(target_region,1))'); 
    target_mask = (x-target_centre(1)).^2 + (y-target_centre(2)).^2 < max_radius^2;
    
    %Compute the orientation hist for the target region
    [tout rout target_orientation_hist] = ...
        weighted_complex_rose(target_orientation_map(target_mask), n_bins);
    
%     subplot(1,2,2); polar(tout, rout); title('Target orientation');
    
    %Compute the correlation between the histograms for all shifts
    orientation_corr = zeros(n_bins,1);
    for o = 1:n_bins
        orientation_corr(o) = ...
            corr2(target_orientation_hist, circshift(template_orientation_hist,o));
    end
    
    %Take rotation angle as the shift that produces maximum correlation
    [dummy rotation_angle] = max(orientation_corr);
    clear dummy;
    display(rotation_angle);
    
    %Convert rotation_angle to radians
    rotation_angle = 2*pi*rotation_angle/n_bins;

end

%Modify the target region using the template
modify_args.MassDilate = 100;
modify_args.Mass = mass_template;
modify_args.TargetImage = target_region;
modify_args.Location = target_centre;
modify_args.Rotate = rotation_angle;
[modified_region target_centre] = mb_transfer_mass_dual_tree(modify_args);
clear modify_args;

% %Now sample a new mass in the modified target region
% sample_args.TargetRegion = modified_region;
% sample_args.TargetCentre = target_centre;
% sample_args.ConditionMode = mass_idx;
% sample_args.Rotation = mass_model.rotations(:,:,mass_idx)*...
%     [cos(rotation_angle) sin(rotation_angle); -sin(rotation_angle) cos(rotation_angle)];
% sample_args.Origin = mass_model.origins(mass_idx);
% sample_args.MassModel = mass_model;
% sample_args.NumModes = 5;
% [mass] = sample_new_mass_in_region(sample_args);
% mass.template_idx = mass_idx;
mass = modified_region;
