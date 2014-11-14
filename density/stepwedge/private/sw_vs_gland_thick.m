% Stepwedge calibration data
% Use this to generate a relationship between gland thickness, breast
% thickness and stepwedge thickness

function [gland_thickness gland_thicknesses] = sw_vs_gland_thick(sw_lookup, breast_thickness, calibration_data, debug_mode)

if nargin < 4
    debug_mode = 0;
end

% calibration is nx3 array where each row contains the breast thickness,
% gland thickness and step height for a single calibration image (this
% shouldn't ever happen from the way calibration is formed, but just in
% case...
[dummy unique_data] = unique(calibration_data(:,1:2), 'rows');
calibration_data = calibration_data(unique_data,:);

%Find unique breast thicknesses in calibration_data
[breast_thicknesses, dummy, bt_idx] = unique(calibration_data(:,1));

n_thick = length(breast_thicknesses);

%Create structure to store gland thicknesses in
step_heights = 0.5:0.5:14;
n_steps = length(step_heights);
gland_thicknesses = zeros(n_thick, n_steps);

%For each breast thickness get the gland thickness and step height pairs from
%calibration data then interpolate to find gland thickness at all step
%heights in the wedge
for ii = 1:n_thick
    
    %Get gland thicknesses and step heights for this breast thickness
    gt = calibration_data(bt_idx == ii, 2)';
    sh = calibration_data(bt_idx == ii, 3)';
    
    %Do 1-D interp to get gland_thickness at each step height in the wedge
    gland_thicknesses(ii,:) = ...
        interp1(sh, gt, step_heights,'linear','extrap');
    
    %If we've extrapolated to any gland thickness larger than the breast
    %thickness this must be wrong, so limit at breast thickness 
    gland_thicknesses(ii, gland_thicknesses(ii,:) > breast_thicknesses(ii))...
        = breast_thicknesses(ii);
end

%We also can't have any negative gland thicknesses!
gland_thicknesses(gland_thicknesses < 0) = 0;

if ~isempty(sw_lookup) && ~isempty(breast_thickness)
    gland_thickness = interp2(...
        repmat(step_heights, n_thick, 1),... %Grid of all step heights
        repmat(breast_thicknesses(:), 1, n_steps),... % 
        gland_thicknesses,...
        sw_lookup, breast_thickness, 'linear');
else
    gland_thickness = [];
end

if debug_mode
    [sh bt] = meshgrid(linspace(0,14,100), linspace(0,max(breast_thicknesses(:)),100));
    
    gt = interp2(...
    repmat(step_heights, n_thick, 1),... %Grid of all step heights
    repmat(breast_thicknesses(:), 1, n_steps),... % 
    gland_thicknesses,...
    sh, bt, 'linear');
    figure('Name', 'Calibration surface');
    mesh(sh,bt,gt);
    xlabel('Step heights (mm)');
    ylabel('Breast thickness (mm)');
    zlabel('Gland thickness (mm)');
    title('Surface fit calibrating gland thickness to step height and breast thickness');
end
    

