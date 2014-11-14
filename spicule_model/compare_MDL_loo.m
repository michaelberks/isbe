% COMPARE_MDL_LOO put mass outlines into shape matrix 
%    [shapes] = compare_MDL_loo(model_shapes, mass_path, npt)
%
%    inputs:
%       mass_files  - Structure listing file names of input masses
%                       each DxN matrix corresponds to a particular shape
%       mass_path   - File path to mass folder
%       n_pts       - number of points to use per shape          
%                       D is the dimensionality of the data points 
%
%    outputs:
%       
%
%    notes:
%
% Copyright: (C) 2006-2008 University of Manchester
% Author: Michael Berks

function [er_pro er_mdl shapes_m shapes_t shapes_rg] =...
    compare_MDL_loo(model_shapes, test_shapes, mass_path, ssv)



    %
    % Model building phase
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Load model shapes and standardise n_pts
    N_model = length(model_shapes);
    [shapes_m.unaligned, mass_areas] =...
        get_shapes_from_masses(mass_files, mass_path, n_pts);
    ss = get_shapes_from_masses(mass_files, mass_path, n_pts, 'mdl');
    
    
    % Make shape model by procrustes alignment (i.e. original method)
    %%%
    [mean_pro, P_pro, a, a,...
    a, a, a, a, a, a, shapes_m.pro] ...
        = shape_model(shapes_m.unaligned, mean(mass_areas), 0.98);
    clear a;
    
    % Align shapes using the MDL model
    %%
    [r_pts] = amb_automodelbuild (ss, 'circle', 'nIterations', 10, 'saveFrequency', 10,...
        'Quiet', 0, 'optimisePose', 0, 'optimiseOrigin', 0, 'nExamplesToOptimise', N_model,...
        'initialAlign', 1, 'initialOriginOptimisation', 1);
    
    % Reshape the reparametrised points
    shapes_m.mdl = zeros(50, 1000);
    for ii = 1:50; shapes_m.mdl(ii,:) = [r_pts(1,:,ii) r_pts(2,:,ii)]; end
    
    % Compute model parameters for MDL shapes
    [mean_mdl, P_mdl] = pca(shapes_m.mdl, size(P_pro,2));
    
    %
    % Testing phase
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Load test shapes and standardise n_pts
    N_test = length(test_shapes);
    shapes_t.unaligned = zeros(N_test, 2*ssv);
    shapes_t.pro = zeros(N_test, 2*ssv);
    shapes_t.mdl = zeros(N_test, 2*ssv);
    for ii = 1:N_test
        temp = load([mass_path, test_shapes(ii).name]);
        mass = temp.mass; clear temp;

        idx = round(linspace(1, length(mass.mass_outline(:,1)), ssv+1));
        idx(end) = []; %ensures first point is not equal to the last point!
        shape_vec = mass.mass_outline(idx,:);
        shapes_t.unaligned(ii,:) = [shape_vec(:,1)', shape_vec(:,2)'];
                
        [dd Z] = procrustes([mean_pro(1:ssv)' mean_pro(ssv+1:end)'],...
            shape_vec);
        shapes_t.pro(ii, :) = [Z(:,1)' Z(:,2)'];
        
        [dd Z] = procrustes([shapes_m.mdl(1,1:ssv)' shapes_m.mdl(1,ssv+1:end)'],...
            shape_vec);
        shapes_t.mdl(ii,:) = [Z(:,1)' Z(:,2)'];
        clear mass shape_vec idx;
    
    
    end
    
    % Compute shape error for two models
    B_pro = P_pro' * (shapes_t.pro - repmat(mean_pro, N_test, 1))';
    B_mdl = P_mdl' * (shapes_t.mdl - repmat(mean_mdl, N_test, 1))';
    
    shapes_rg.pro = repmat(mean_pro, N_test, 1) + (P_pro*B_pro)';
    shapes_rg.mdl = repmat(mean_mdl, N_test, 1) + (P_mdl*B_mdl)';
    
    XX = shapes_t.unaligned(:,1:ssv); YY = shapes_t.unaligned(:,ssv+1:end);
    Xp = shapes_rg.pro(:,1:ssv); Yp = shapes_rg.pro(:,ssv+1:end);
    Xm = shapes_rg.mdl(:,1:ssv); Ym = shapes_rg.mdl(:,ssv+1:end);
    
    er_pro = mean(sqrt((XX - Xp).^2 + (YY- Yp).^2), 2);
    er_mdl = mean(sqrt((XX - Xm).^2 + (YY- Ym).^2), 2);
    
    
    
    