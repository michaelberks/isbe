function [mass_spline] = subtract_mass_dual_tree(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% author:   Michael Berks
% date:     11/04/2006  11:50
%
% function: Use TPS to interpolate grey-levels from landmark pts in 
% ROI surrounding mass and subtract mass background. Ouput structure of ROIs
% of subtracted mass gray levels
%
%inputs:
%   args.MassList: list of mass files to apply background estimation to
%   args.Path: directory name where masses are located
%   args.n1: no. of pixels from shape border to start of landmark pts
%   args.n2: no. of pixels from shape border to end of landmark pts
%   args.Spacing: no. of pixels between landmark pts in x and y directions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Unpack the arguments:
args = u_packargs(varargin, '0', ... % non strict mode
		  {...  % The mandatory arguments
          'MassList'...
          },... % The optional arguments
          'MassDir', 'C:\isbe\dev\annotations\',...
          'Level', 5,...
          'n1', 50,...
          'n2', 50, ...
          'Spacing', 1, ...
          'Sigma', 11, ...
          'Green', 'biharmTPS', ...
          'FilterMethod', 0, ...
          'Plot', 0 ...
          );
      
%Get number of masses 
N = length(args.MassList);

%reduce size of n1 and n2 to match level scaling
% args.n1 = args.n1 / 2^(args.Level - 1);
% args.n2 = args.n2 / 2^(args.Level - 1);

for ii = 1:N
    mass = u_load([args.MassDir, args.MassList(ii).name]);
    
    %Compute DT-CWT decomp of mass region, and set region at the
    %args.Level-th scaling coefficents
    [mass_dt mass_sc] = dtwavexfm2(mass.mass_ROI, args.Level);
    mass_ROI = mass_sc{args.Level};
    mass_ROI_fine = mass_sc{args.Level-1};
    
    
    
    %extract shape border input structure and reduce size to match level of
    %scaling coefficients
    shape_vec = mass.mass_outline / 2^(args.Level - 1);
    shape_vec_fine = mass.mass_outline / 2^(args.Level - 2);
    
    shape_bw = roipoly(mass_ROI, shape_vec(:,1),...
        shape_vec(:,2));
    
    shape_bw_fine = roipoly(mass_ROI_fine, shape_vec_fine(:,1),...
        shape_vec_fine(:,2));
    
    for jj = 1:(2*args.n1 - 1)
        shape_bw_fine = imdilate(shape_bw_fine, strel('disk', 1));
    end   
    
    display('start iterating');
    under_ssd = 3;
    while under_ssd > 1
        switch args.FilterMethod
            case 0
                smooth_ROI = mass_ROI;
            case 1
                smooth_ROI = imfilter(mass_ROI, gauss_filt, 'symmetric');
                f_string = 'Gaussian';
            case 2
                smooth_ROI = medfilt2(mass_ROI, [args.Sigma args.Sigma], 'symmetric');
                f_string = 'Median';
            case 3
                smooth_ROI = wiener2(mass_ROI, [args.Sigma args.Sigma], 'symmetric');
                f_string = 'Wiener';
        end

        [bg_estimates p_list] = spline_estimation(smooth_ROI, shape_bw,...
            args.n1, args.n2, args.Spacing, args.Green);
    
        %Subtract interpolated values from background and create new ROI of
        %subtracted gray-levels (zero outside of shape)
        
        bg_diff = double(mass_ROI(p_list)) - bg_estimates';
        total_ssd = sum(bg_diff.^2) / length(p_list);
        under_ssd = sum(bg_diff(bg_diff > 0).^2) / length(p_list);
        display(['total = ', num2str(total_ssd), ' under = ', num2str(under_ssd)]);
        
        if under_ssd > 0.5*total_ssd   
            mass_ROI(p_list) = uint8(bg_estimates);
        else
            break;
        end
    end
    display('end iterating');
    display('');
    
    %Reconstruct the mass region from the dual-tree
    mass_dt{args.Level+1} = mass_ROI;
    clear mass_ROI bg_estimates;
    
    %Now the clever bit - after reconstructing the coarsest level, swap
    %back the original scaling coefficients in the region outside the mass,
    %so that the global change is reduced
    [mass_spline] = dtwaveifm2(mass_dt, args.Level);
    mass_spline(~shape_bw_fine) = mass_ROI_fine(~shape_bw_fine);
    
    mass_dt(args.Level+1) = [];
    mass_dt{args.Level} = mass_spline;
    
    %Now reconstruct the full tree;
    [mass_spline] = dtwaveifm2(mass_dt);
    
    
    
    if args.Plot
        figure('Name', ['Mass', zerostr(ii,3)]);
        imagesc(mass_spline); axis image; %colormap gray;
        hold on;
        %plot(l_c, l_r, 'rx');        
    end
    
    save(['C:\isbe\dev\misc\anno_bg\', args.MassList(ii).name], 'mass_spline');
        
end

% display(['Finished subtracting background:', ...
%                 ' args.n1 = ', num2str(args.n1), ...
%                 ' args.n2 = ', num2str(args.n2), ...
%                 ' args.Spacing = ', num2str(args.Spacing), ...
%                 ' args.Sigma = ', num2str(args.Sigma),...
%                 ' filter = ', f_string]);