% SAMPLE_NEW_MASSES
%    [] = sample_new_masses(varargin)
%
% SAMPLE_NEW_MASSES generates and saves a set of new masses from a mass appearance
% model
%
% SAMPLE_NEW_MASSES uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%   'MassModel' - Mass appearance model
%   'NumberOfMasses' - Number of masses to generate 
%   'SavePath' - path to directory in which to save masses
%
% Optional Arguments:
%   'SaveName' - file name (to be appended with idx) of mass. Deafult: 'new_mass'
%   'Plot' - {0,1} display each new mass in figure
%
% Return Value:
%
% Copyright:(C) 2006-2008 University of Manchester
% Author:   Michael Berks
% Date:     08/06/2006  09:52
function [mass] = sample_scale_mass_in_region(varargin)

args = u_packargs(varargin, '0', ... % strict mode
		  {'MassModel',...
          'TargetRegion',...
          'TargetCentre',...
          'ScaleParameter',...
          },...
          'ConditionMode', [],...
          'ConditionLevel', 1,...
          'SavePath', [],...
          'SaveName', 'new_mass',... % The optional arguments
          'Plot', 0);
      
%
% Get args from args.MassModel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mean_shape  = args.MassModel.mean_shape;
P_shape     = args.MassModel.P_shape;
% B_shape     = args.MassModel.B_shape;
L_shape     = args.MassModel.L_shape;
mean_scale  = args.MassModel.mean_scale;
P_scale     = args.MassModel.P_scale;
% B_scale     = args.MassModel.B_scale;
% L_scale     = args.MassModel.L_scale;
mean_tex    = args.MassModel.mean_tex;
P_tex       = args.MassModel.P_tex;
% B_tex       = args.MassModel.B_tex;
L_tex       = args.MassModel.L_tex;


% mean_com      = args.MassModel.mean_com;
P_com         = args.MassModel.P_com;
% B_com         = args.MassModel.B_com;
L_com         = args.MassModel.L_com;

W_shape     = args.MassModel.W_shape;
W_tex       = args.MassModel.W_tex;
W_scale     = args.MassModel.W_scale;
% W_n         = args.MassModel.W_n;

C_nm        = args.MassModel.C_nm;
C_mm_inv    = args.MassModel.C_mm_inv;

mean_shape_pl = args.MassModel.mean_shape_pl;
size_shape_vec = length(mean_shape) / 2;

k_shape     = length(L_shape);
k_tex       = length(L_tex);

scaled_mean_shape = args.MassModel.scaled_mean_shape;
  
%Define source points for TPS - as row vectors
s_x = scaled_mean_shape(1:size_shape_vec);% + mean_centre(1);
s_y = scaled_mean_shape(size_shape_vec+1:end);% + mean_centre(2);

%Define points to be interpolated by TPS - as row vectors
i_x = mean_shape_pl(:,1)';
i_y = mean_shape_pl(:,2)';
        
%tps_L_inv = tps_weights(s_x, s_y);

%
% Generate specified number of masses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Make the output directory if it doesn't already exist
if 0
if ~isdir(args.SavePath)
    mkdir(args.SavePath);
end

% ensure the directory names are filesep-terminated
if ~strcmp(args.SavePath(end), filesep)
	args.SavePath = [args.SavePath filesep];
end
end

% sample new combined appearance vectors - assume normal distribution
% of modes
%%%

% Compute new combined shape and texture vector, conditioned on scale

mu = -C_mm_inv*C_nm*[args.ScaleParameter];
[e_vec, e_val] = eig(C_mm_inv);
new_com = e_vec*sqrtm(e_val)*randn(length(L_com), 1) + mu;

%If we've been given the first mass mode, substitute this here:
if ~isempty(args.ConditionMode)
    new_com(1:args.ConditionLevel) = args.MassModel.B_com(1:args.ConditionLevel,args.ConditionMode);
end

Q_shape = P_com(1:k_shape,:); 
Q_tex = P_com(k_shape+1:k_shape + k_tex,:);

new_shape = mean_shape + (P_shape*inv(W_shape)*Q_shape*new_com)';
new_tex = mean_tex + (P_tex*inv(W_tex)*Q_tex*new_com)';
new_scale = mean_scale + P_scale*args.ScaleParameter / W_scale;

new_shape = new_shape/new_scale;
new_shape = reshape(new_shape, [], 2);

new_shape(:,1) = new_shape(:,1) + args.TargetCentre(1);
new_shape(:,2) = new_shape(:,2) + args.TargetCentre(2);

new_bw = roipoly(args.TargetRegion, new_shape(:,1), new_shape(:,2));

%     for jj = 1:49; new_bw = imdilate(new_bw, strel('disk', 1)); end
[new_label, num_obj] = bwlabel(new_bw, 4);    

if num_obj == 1
    rp = regionprops(new_label, 'PixelList'); clear new_label;
    new_shape_pl = rp.PixelList; clear rp;
    clear rp new_bw;

    % Compute TPS warp to map from mean to new shape
    %%%%

    %Define displacement to target points
    z_x = new_shape(:,1)';
    z_y = new_shape(:,2)';

    %Compute displacement of interpolated points        
    T = geom_alignpoints_2d([z_x; z_y], [s_x; s_y], [],...
        'transform', 'spline');
    [pts] = geom_transformpoints([i_x; i_y], T);
    f_x = pts(1,:);
    f_y = pts(2,:);

    display('completed warping');
    
    % Create new shape pixel list
    %%%%%%%%%%

    new_shape_tex = griddata(f_x, f_y, new_tex,...
        new_shape_pl(:,1), new_shape_pl(:,2));
    new_shape_tex(isnan(new_shape_tex)) = 0;
    new_shape_ROI = zeros(size(args.TargetRegion));
    new_shape_ROI(sub2ind(size(args.TargetRegion), new_shape_pl(:,2), new_shape_pl(:,1)))...
        = new_shape_tex;

    display('completed grid');

    

    mass.mass_outline = new_shape;
    mass.subtract_ROI = new_shape_ROI;
    mass.mass_ROI = args.TargetRegion + new_shape_ROI;
    mass.mass_centroid = args.TargetCentre;
    
    if ~isempty(args.SavePath)
        mass_name = [args.SavePath, args.SaveName];
        save(mass_name, 'mass');
    end
    if args.Plot

        %figure('WindowStyle', 'docked');
        %imagesc(new_shape_ROI); axis image; colormap(gray(256));
        %hold on
        %plot(new_shape(:,1),new_shape(:,2), 'y','LineWidth', .5);
        %plot(f_x(1:3:end), f_y(1:3:end), 'rx', 'MarkerSize', .5);
        figure; 
        subplot(1,3,1); imagesc(uint8(mass.subtract_ROI)); axis image; colormap(gray(256));
        subplot(1,3,2); imagesc(uint8(mass.mass_ROI)); axis image; colormap(gray(256));
        subplot(1,3,3); imagesc(uint8(mass.mass_ROI-mass.subtract_ROI)); axis image; colormap(gray(256));
    end
    clear new_shape_tex new_shape_pl new_shape_ROI new_tex;
    
else
    display('Bad shape created');
end
