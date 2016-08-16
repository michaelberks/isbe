function [res, Jrows, Jcols, Hrows, Hcols] = ...
    smoothness_error_vec(var_vec, sz_vec, observations, ...
                         f_clear_persistent)
                     
if (nargin==0 && nargout==0), test(); return; end
          
persistent Us;
persistent Vs;
persistent Ts;

if ~exist('f_clear_persistent','var'), f_clear_persistent = false; end

% Reset the persistent variables and do nothing more.
if f_clear_persistent
    clear('Us'); 
    return; 
end

variables = ncm_unpack(var_vec, sz_vec);
[p, wts, stackSize, uu, vv] = get_params_from_vec(observations, variables);

res = [];

%% Data smoothness
if isempty(Us)
    m = stackSize(1);
    n = stackSize(2);
    
    tau_hv = stackSize(3) / 6;
    tau_d = 0;

    if ~isempty(observations.presence)
        img = 0 * double(observations.presence > 0.5);
    else
        % Since spatial things are tied to the first frame, use this instead of
        % the mean (unless you want to register them all)
        img = 0 * mean(observations.imgStack, 3);
    end

    [Us, Vs, Ts] = create_spatial([m,n], 1:n, 1:m, ...
                                  tau_hv, tau_d, img);
end

res = [res;
       Us*uu(:) + Ts;
       Vs*vv(:) + Ts];

   
% %% Displacement smoothness
% p = reshape_variables(p);
% diffs = p.displacements(2:end,:) - p.displacements(1:end,:);
% diffs_squared = diffs.^2;
% 
% res = [res;
%        diffs_squared(:)]; 

% Compute locations of nonzero elements of the Jacobian adn Hessian 
% for this vector of residuals
if (nargout >= 3)
    [Jrows, Jcols] = jacobian_elements(var_vec, sz_vec, Us);
end
if (nargout >= 5)
    [Hrows, Hcols] = hessian_elements(var_vec, sz_vec, Us, Jrows, Jcols);
end

% return


function [rows, cols] = jacobian_elements(var_vec, sz_vec, Us)

[s_cols, s_rows] = find(Us');
nConsSpatial = size(Us,1);

[vs_names, vt_names] = varnames();
vnames = [vs_names vt_names];

rows = [];
cols = [];

used_var_inds = find(sz_vec ~= 0)';
for v = used_var_inds
    varname = vnames{v};
    
    switch varname
        case {'uu', 'vv'},
            % Variable that affect *either* the smoothness in U *or* the
            % smoothness in V but not both.
            
            var_ind = find(strcmp(varname, vnames(sz_vec ~= 0)));

            if strcmp(varname,'uu')
                row_offset = 0;
            else
                row_offset = nConsSpatial;
            end
                
            rows = [rows(:); row_offset + s_rows(:)];

%             cols = (col_spacing * (s_cols(:)-1)) + iv;
%                  = (col_spacing * s_cols(:)) + (iv - col_spacing);
            col_spacing = sum(sz_vec(1:length(vs_names)) ~= 0);
            col_offset = var_ind - col_spacing;
            cols = [cols(:); col_offset + (col_spacing*s_cols(:))];

        case {'Q', 'F', 'theta'},
            % Variables that affect *both* the smoothness in U *and* the
            % smoothness in V.
            
            var_ind = find(strcmp(varname, vnames(sz_vec ~= 0)));

            col_spacing = sum(sz_vec(1:length(vs_names)) ~= 0);
            
            for row_offset = [0 nConsSpatial]
                rows = [rows(:); row_offset + s_rows(:)];
                
%                 cols = (col_spacing * (s_cols(:)-1)) + iv;
%                      = (col_spacing * s_cols(:)) + (iv - col_spacing);
                col_offset = var_ind - col_spacing;
                cols = [cols(:); col_offset + (col_spacing*s_cols(:))];
            end
            
        case {'f'},
            % Temporal flow variations do not affect the smoothness (much).
            % (They do in the sense that we should technically scale the
            % penalties accordingly but I'm ignoring that for now.)
            
        case {'displacements'}
            % Displacements really don't affect the smoothness of the flow
            % field.
    end
end


function [rows, cols] = hessian_elements(var_vec, sz_vec, Us, Jrows, Jcols)

nv = sum(sz_vec);
nr = 2*size(Us,1);

J = sparse(Jrows,Jcols, 1, nr,nv);

[rows, cols] = find(J'*J);

return

% Find rows corresponding to spatial penalties
[s_cols, s_rows] = find(Us');

vnames = varnames();

%% Fill in the Hessian pattern

% i.e., an [nv x nv] matrix, H, where H(i,j)=1 indicates that at
% least one of the observations influences both variable i and variable j.

rows = [];
cols = [];

u_found = false;
for v = find(sz_vec ~= 0)
    varname = vnames{v};
    
    switch varname
        case ('Q'),
            % Do nothing for now.
            
        case ('u'),
            u_found = true;
            
        case ('v'),
            if (u_found)
                % Both u and v at a given pixel are influenced via the data
                % (brightness) error

                % Reshape into a matrix whose columns correspond to pairs of row/column
                % indices. Then repeat appropriately so that the indices give a
                % block-diagonal matrix of 2x2 submatrices.
                tmp = reshape(1:sum(sz_vec(1:nvs)), 2,[]);
                rows = [rows; kron(tmp, [1;1])];
                cols = [cols; kron([1;1], tmp)];

                % Two u's and v's at neighbouring pixels are influenced by the spatial
                % smoothness error

                % s_cols are ordered such that every consecutive pair correspond to the
                % same observation
                v1 = s_cols(1:2:end);
                v2 = s_cols(2:2:end);

                % U and V lie on interleaved columns
                % Matrix is symmetric so any additions to (row,col) also apply to
                % (col,row) elements.
                rows = [rows(:); 2*v1(:)-1; 2*v1(:); 2*v2(:)-1; 2*v2(:)];
                cols = [cols(:); 2*v2(:)-1; 2*v2(:); 2*v1(:)-1; 2*v1(:)];
            end
            
        case {'f'},
            % f influences every pixel in the same but nothing more
            rng = 1:nImages;
            rng = rng(:) + sum(sz_vec(1:nvs));
            rows = [rows; rng(:)];
            cols = [cols; rng(:)];
            
        case {'displacements'}
            % displacements two frames apart affect every pixel in the
            % intermediate frame.

            % Diagonal plus two off-diagonals on either side
            r = [1:nImages-1 1:nImages-2 1:nImages-3 2:nImages-1 3:nImages-1];
            c = [1:nImages-1 2:nImages-1 3:nImages-1 1:nImages-2 1:nImages-3];
           
            off2 = sum(sz_vec(1:nvs+1));
            
            if (0)
                % FIXME: This needs completing correctly before we can
                % optimize over f and displacements together.
                
                % f (at time t) and displacement (at times t-1 and t+1) are
                % both influenced by the pixels in frame t
                off1 = sum(sz_vec(1:nvs));
                % Add horizontal blocks
                rows = [rows; r(:)+off1; r(:)+off1]; 
                cols = [cols; c(:)+off2; c(:)+off2+nImages-1];
                % Add vertical blocks
                rows = [rows; r(:)+off2; r(:)+off2+nImages-1];
                cols = [cols; c(:)+off1; c(:)+off1];
            end

            % Offset to first column after variables.f
            r = r(:) + off2;
            c = c(:) + off2;

            rows = [rows; r; r; r+nImages-1; r+nImages-1];
            cols = [cols; c; c+nImages-1; c; c+nImages-1];
    end
end

% Temporal variable dependencies
% Every pair of temporal variables is influenced by the same pixels so
% create a block of ones.
[c,r] = meshgrid(sum(sz_vec(1:nvs))+1:sum(sz_vec), 1:sum(sz_vec(1:nvs)));
rows = [rows; r(:); c(:)];
cols = [cols; c(:); r(:)];


%% Test script
function test()
clc;

[variables, observations, var_vec, sz_vec] = create_dummy_variables(false);

% Clear the persistents first
smoothness_error_vec([],[],[], true);

[res, Jrows, Jcols] = smoothness_error_vec(var_vec, sz_vec, observations);

nv = sum(sz_vec);
nr = length(res);

JP = sparse(Jrows,Jcols, 1, nr,nv);
HP = JP'*JP;

figure(1); clf;
    spy(JP);
      
figure(2); clf;
    spy(HP);
