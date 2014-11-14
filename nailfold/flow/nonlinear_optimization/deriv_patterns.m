%% Generate an image of the Jacobian, dz/dx, using known dependencies.
function [J, H] = deriv_patterns(var_vec, sz_vec, imsz, ...
                                 observationMask, f_clear_persistent)
if (nargin==0 && nargout==0), test(); return; end
                               
if ~exist('observationMask','var'), observationMask = []; end
if ~exist('f_clear_persistent','var'), f_clear_persistent = false; end
                               
persistent s_rows;
persistent s_cols;
persistent nConsSpatial;

if f_clear_persistent
    clear('s_rows','s_cols','nConsSpatial');
    return; 
end

if isempty(s_rows)
    m = imsz(1);
    n = imsz(2);
    Us = create_spatial([m,n], 1:n, 1:m, ...
                        1/6, 0);
                    
    % Find rows corresponding to spatial penalties
    [s_cols,s_rows] = find(Us');
    nConsSpatial = size(Us,1);
end

nPixelsPerImage = imsz(1)*imsz(2);
nImages = imsz(3);

nv = length(var_vec);
v_nz = (sz_vec ~= 0);

[vs, vt] = varnames();

nvs = length(vs);
vs_nz = v_nz(1:nvs);
nvs_nz = sum(vs_nz);

nvt = length(vt);
vt_nz = v_nz(nvs+1:nvs+nvt);
nvt_nz = sum(vt_nz);

%% Fill in the Hessian pattern

% i.e., an [nv x nv] matrix, H, where H(i,j)=1 indicates that at
% least one of the observations influences both variable i and variable j.

rows = [];
cols = [];

% Spatial variable dependencies.
if (nvs_nz > 0)
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

% Temporal variable dependencies
% Every pair of temporal variables is influenced by the same pixels so
% create a block of ones.
[c,r] = meshgrid(sum(sz_vec(1:nvs))+1:sum(sz_vec), 1:sum(sz_vec(1:nvs)));
rows = [rows; r(:); c(:)];
cols = [cols; c(:); r(:)];

for v = find(vt_nz)'
    varname = vt{v};
    switch varname
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

H = sparse(rows,cols, true, nv,nv);


%% Fill in the Jacobian pattern
if isempty(observationMask)
    nObservations = nImages*nPixelsPerImage;
else
    obsInds = find(observationMask);
    nObservations = length(obsInds);
end

rows = [];
cols = [];

if (nvs_nz > 0)
    % Find offsets for spatial smoothness penalties.
    i = find(strcmp(vs, 'uu'));
    uu_offset = sum(sz_vec(1:i-1));

    if isempty(observationMask)
        nObservations = nImages*nPixelsPerImage;

        rowmat = zeros(nImages, nPixelsPerImage*nvs_nz);
        rowmat(1,:) = kron(1:nPixelsPerImage, ones(1,nvs_nz));
        for r = 2:nImages
            rowmat(r,:) = rowmat(r-1,:) + nPixelsPerImage;
        end

        colmat = repmat(1:nPixelsPerImage*nvs_nz, [nImages,1]);
        
        rows = [rowmat(:); 
                nObservations + s_rows(:);
                nObservations + nConsSpatial + s_rows(:)];

        % Columns are interleaved for better arrowhead structure
        cols = [colmat(:); 
                uu_offset + (2*s_cols(:))-1;
                uu_offset + (2*s_cols(:))];
    else
        obsInds = find(observationMask);
        nObservations = length(obsInds);
        
        rowmat = zeros(size(observationMask));
        rowmat(obsInds) = 1:nObservations;
        rowmat = reshape(rowmat, [nPixelsPerImage, nImages]);
        rowmat = repmat(rowmat', [nvs_nz, 1]);
        
        nz_inds = find(rowmat ~= 0);

        % Data terms
        rows = rowmat(nz_inds);
        cols = ceil(nz_inds / nImages);

        % Spatial terms
        % (Columns are interleaved for better arrowhead structure)
        
        % Smoothness in U
        rows = [rows; nObservations + s_rows(:);];
        cols = [cols; uu_offset + (2*s_cols(:))-1];

        % Smoothness in V
        rows = [rows; nObservations + nConsSpatial + s_rows(:)];
        cols = [cols; uu_offset + (2*s_cols(:))];
    end
end

% Temporal variable dependencies.
col = nPixelsPerImage*nvs_nz + 1;
for v = find(vt_nz)'
    varname = vt{v};
    switch varname
        case {'f'},
            % f affects only those frames to which it applies directly
            if isempty(observationMask)
            else
                rowmat = zeros(size(observationMask));
                rowmat(obsInds) = 1:nObservations;
                
                for f = 1:nImages
                    row_submat = rowmat(:,:,f);
                    nz_inds = (row_submat ~= 0);
                    rows = [rows; row_submat(nz_inds)];
                    cols = [cols; col*ones(sum(nz_inds(:)),1)];
                    col = col + 1;
                end
            end
            
        case {'displacements'}
            % displacements affect frames immediately before and after 
            % (via the temporal gradients)
            
            if isempty(observationMask)
            else
                rowmat = zeros(size(observationMask));
                rowmat(obsInds) = 1:nObservations;
                
                % Displacement doesn't exist for frame 1...
                % ...but for the final frame, it affects both nImages-1 and
                % nImages
                for f = 2:nImages
                    frames = min(max(f-1:f+1, 1), nImages);
                    frames = unique(frames);
                    
                    row_submat = rowmat(:,:,frames);
                    nz_inds = (row_submat ~= 0);

                    % Displacement in both X and Y
                    % (Remember that columns are interleaved)
                    rows = [rows; row_submat(nz_inds);
                                  row_submat(nz_inds)];
                    cols = [cols; col*ones(sum(nz_inds(:)),1);
                                  (col+nImages-1)*ones(sum(nz_inds(:)),1)];
                    col = col + 1;
                end
            end
    end
end

% Define total number of residuals
nr = nObservations + ...  % Data penalties
     2*nConsSpatial;      % Smoothness penalties

J = sparse(rows,cols, true, nr,nv);


% % Sad face...
% spatial_var_inds = sum(sz_vec(1:nvs))+1:size(H,2);
% H(spatial_var_inds,spatial_var_inds) = 1;
% J(:,spatial_var_inds) = 1;


function test()
clc;

[variables, observations, var_vec, sz_vec] = create_dummy_variables();
imsz = size(observations.imgStack);

deriv_patterns([],[],[],[], true);
[JP, HP] = deriv_patterns(var_vec, sz_vec, imsz, observations.obsMask);

figure(1); clf;
    subplot(2,1,1); spy(HP);
    subplot(2,1,2); spy(JP'*JP);
      
figure(2); clf;
    spy(JP);

    