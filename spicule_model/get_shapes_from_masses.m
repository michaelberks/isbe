% GET_SHAPES_FROM_MASSES put mass outlines into shape matrix 
%    [shapes] = get_shapes_from_masses(mass_files, ...)
%
%    inputs:
%       mass_files  - Structure listing file names of input masses
%                       each DxN matrix corresponds to a particular shape
%       mass_path   - File path to mass folder
%       n_pts       - number of points to use per shape          
%                       D is the dimensionality of the data points
%       type        - {'standard'} - shape matrix is (N x 2*n_pts) 
%                      'mdl' - shape matrix is (2 x n_pts x N) 
%
%    outputs:
%       shapes      - the matrix of shapes
%
%    notes:
%
% Copyright: (C) 2006-2008 University of Manchester
% Author: Michael Berks

function [shapes, mass_areas] =...
    get_shapes_from_masses(mass_files, mass_path, n_pts, type)


N = length(mass_files);

switch type
    case 'standard'
        shapes = zeros(N, 2*n_pts);
        mass_areas = zeros(N,1);
        for ii = 1:N
            %load mass
            temp = load([mass_path, mass_files(ii).name]);
            mass = temp.mass; clear temp;
            
            idx = round(linspace(1, size(mass.mass_outline, 1), n_pts+1));
            %ensure first point is not equal to the last point
            shapes(ii,:) = [mass.mass_outline(idx(1:end-1),1)',...
                            mass.mass_outline(idx(1:end-1),2)'];
            mass_areas(ii) = mass.mass_area;
            clear mass idx;
        end

    case 'mdl'
        shapes = zeros(2, n_pts, N);
        for ii = 1:N
            %load mass
            temp = load([mass_path, mass_files(ii).name]);
            mass = temp.mass; clear temp;

            idx = round(linspace(1, size(mass.mass_outline, 1), n_pts+1));
            %ensure first point is not equal to the last point
            shapes(:,:,ii) = [mass.mass_outline(idx(1:end-1),1)';...
                            mass.mass_outline(idx(1:end-1),2)'];
            clear mass idx;
        end
end
