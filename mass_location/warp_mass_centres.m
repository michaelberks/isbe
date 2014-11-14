function [mass_centres] = warp_mass_centres(mean_shape, mass_list, mlo, debug_mode)
%WARP_MASS_CENTRES *Insert a one line summary here*
%   [mass_centres] = warp_mass_centres(mean_shape,mass_list,debug_mode)
%
% Inputs:
%      mean_shape- *Insert description of input variable here*
%
%      mass_list- *Insert description of input variable here*
%
%      debug_mode- *Insert description of input variable here*
%
%
% Outputs:
%      mass_centres- *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 17-Nov-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester 

%Set default parameters
if nargin < 4
    debug_mode = 1;
end

%Get number of masses
n_masses = length(mass_list);

%Pre-allocate space for warped mass centres
mass_centres = zeros(n_masses, 2);

%Annotation directory
anno_dir = 'C:\isbe\dev\annotations\';

%Warp each breast shape to the breast mean, and in doing so compute the
%warped value of the mass centroid
for ii = 1:n_masses
    
    %Load in annotated mass structure
    anno = u_load([anno_dir mass_list(ii).name]);
    
    %Load in mammogram
    mam = imread(anno.name);
    
    %Get mass border
    mass_border = [anno.C1+anno.mass_outline(:,1)-1 anno.R1+anno.mass_outline(:,2)-1];
    
    %Build segmentation from mass name and load segmentation
    seg_list(1).name = ['o', mass_list(ii).name(3:11), '_segmentation.mat'];
    segmentation = u_load(...
        ['C:\isbe\dev\segmentation\breast_borders\', seg_list(1).name]);
    
    %If debugging display mammogram
    if debug_mode
        [large_breast_border] = segment_breast_resize(size(mam), segmentation);
        figure; 
        subplot(1,3,1); imagesc(mam); axis image; colormap(gray(256)); hold on;
        plot(mass_border(:,1), mass_border(:,2));
        plot(large_breast_border(:,1), large_breast_border(:,2), 'r');
    end
    
    nr = segmentation.size(1);
    nc = segmentation.size(2);
    
    %Resize mass border to fit in segmentation
    r = size(mam, 1); clear mam;
    resize_factor = nr / r;
    mass_border = (mass_border - 1)*resize_factor + 1;
    
    %Work out mass centroid w.r.t. resized breast border
    [bw] = poly2mask(mass_border(:,1), mass_border(:,2), nr, nc);
    [yy xx] = find(bw);
    mass_centre = [mean(xx) mean(yy)];
    clear xx yy bw
    
    if ~isempty(strfind(mass_list(ii).name, 'R'));
        mass_border(:,1) = nc - mass_border(:,1) + 1;
        mass_centre(:,1) = nc - mass_centre(:,1) + 1;
    end
    
    %Get breast border as a standard CC/MLO shape
    if mlo
        [breast_shape dummy T] = get_mlo_shape(seg_list, 50);
        breast_shape = reshape(breast_shape, [], 2);
        mass_border = mass_border*T.rot - repmat(T.t,size(mass_border,1),1);
        mass_centre = mass_centre*T.rot - repmat(T.t,size(mass_centre,1),1);
    else
        breast_shape = reshape(get_cc_shape(seg_list, 50),[],2);
    end    
    
    %Align breast shape to the target mean shape
    [dd breast_shape_a t] = mb_procrustes(mean_shape, breast_shape);
    
    %Map mass border and mass centre into the space of the aligned shape
    mass_border_a = t.b*mass_border*t.T + repmat(t.c(1,:), size(mass_border,1),1);
    mass_centre_a = t.b*mass_centre*t.T + repmat(t.c(1,:), size(mass_centre,1),1);
    
    %Now thin-plate spline warp each breast shape to the target mean shape
    
    %Breast shape form the src points
    s_x = breast_shape_a(:,1)';
    s_y = breast_shape_a(:,2)';
    
    %Mean shape forms the target displacements
    z_x = mean_shape(:,1)';
    z_y = mean_shape(:,2)';
    
    %Build geometric transform
    T = geom_alignpoints_2d([z_x; z_y], [s_x; s_y], [],...
        'transform', 'spline');    
    
    %Warp the breast shape, mass border and mass centre to the mean shape
    %space
    [mass_centre_w] = geom_transformpoints(mass_centre_a', T)';
    [mass_border_w] = geom_transformpoints(mass_border_a', T)';
    [breast_shape_w] = geom_transformpoints(breast_shape_a', T)';        
    clear T;
    
    mass_centres(ii,:) = mass_centre_w;
    
    %If debugging display aligned and warped breast shape and mass
    if debug_mode
        subplot(1,3,2); axis ij equal; hold on;
        plot(breast_shape(:,1), breast_shape(:,2), 'b');
        plot(mass_border(:,1), mass_border(:,2), 'r');
        plot(mass_centre(1), mass_centre(2), 'r+');
%         plot(breast_shape_a(:,1), breast_shape_a(:,2), 'b');
%         plot(mass_border_a(:,1), mass_border_a(:,2), 'r');
%         plot(mass_centre_a(1), mass_centre_a(2), 'r+');
 
        subplot(1,3,3); axis ij equal; hold on;
        plot(mean_shape(:,1), mean_shape(:,2), 'b');
        plot(breast_shape_w(:,1), breast_shape_w(:,2), 'r:');
        plot(mass_border_w(:,1), mass_border_w(:,2), 'r');
        plot(mass_centre_w(1), mass_centre_w(2), 'r+');
    end
end