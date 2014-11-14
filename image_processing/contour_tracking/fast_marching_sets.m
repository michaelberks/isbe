function [Tr Fr] = fast_marching_sets(I_ori, x1, y1)

max_dist = 100;
[nrows ncols] = size(I_ori);

nb_offset = [-1 0; 1 0; 0 -1; 0 1];

%Pre-aloocate up maps and masks
alive_mask = false(nrows, ncols);
near_mask = false(nrows, ncols);

x_coords = repmat(1:ncols, nrows, 1);
y_coords = repmat( (1:nrows)', 1, ncols);
to_fill_mask = ((x_coords - x1).^2 + (y_coords - y1).^2) < max_dist^2;

Tr = inf(nrows, ncols);
Fr = zeros(nrows, ncols);
Tr(y1, x1) = 0;
Fr(y1, x1) = 1;

xi = x1;
yi = y1;

%Precompute lookup tables for computing normals to the front
bin_scores = 2.^(0:7);
bin_mask = true(3); bin_mask(2,2) = 0;
[d_fx, d_fy, nearest_x, nearest_y] = precompute_normals();

while any(to_fill_mask(:))
    %Add xi, yi to alive mask and remove it from the near map and fill map
    alive_mask(yi, xi) = 1;
    near_mask(yi, xi) = 0;
    to_fill_mask(yi, xi) = 0;
    
    %Loop through the neighbours
    for i_nb = 1:4
        xr = xi + nb_offset(i_nb,1);
        yr = yi + nb_offset(i_nb,2);
        
        %If they're not in the alive mask compute (or recompute F and T)
        if ~alive_mask(yr, xr)
            
            %add them to the near map (it doesn't matter if they're already
            %there)
            near_mask(yr, xr) = 1;
            
            %Get the neighbours
            nb_r = yr-1:yr+1;
            nb_c = xr-1:xr+1;
    
            %Work out which of them are alive
            local_alive = alive_mask(nb_r, nb_c);
            
            %Use the precomputed lookup to get the normal direction and
            %nearest pixel r'
            bin_lookup = sum(bin_scores(local_alive(bin_mask))) + 1;
            nrx = d_fx(bin_lookup);
            nry = d_fy(bin_lookup);
            nearest_xi = nearest_x(bin_lookup);
            nearest_yi = nearest_y(bin_lookup);
            
            xrp = xr + nearest_xi;
            yrp = yr + nearest_yi;
            
            %Compute Fr and Tr
            e1rp = I_ori(yrp, xrp);
            
            Fr(yr, xr) = min(Fr(yrp, xrp), abs(real(e1rp)*nrx - imag(e1rp)*nry));
            Tri = Tr(yrp, xrp) + 1/Fr(yr, xr);
            
            if Tri < Tr(yr, xr)
                Tr(yr, xr) = Tri;
            end
        end
    end
    [min_val, min_idx] = min(Tr(near_mask));
    x_coords_i = x_coords(near_mask);
    y_coords_i = y_coords(near_mask);
    xi = x_coords_i(min_idx);
    yi = y_coords_i(min_idx);
end 


function [d_fx, d_fy, nearest_x, nearest_y] = precompute_normals()

bin_str = zeros(1,8);
bin_mask = true(3); bin_mask(2,2) = 0;
bin_map = false(3);



d_fx = zeros(256,1);
d_fy = zeros(256,1);
nearest_x = zeros(256,1);
nearest_y = zeros(256,1);

cx = [-1 0 1; -1 0 1; -1 0 1];
cy = cx';

for ii = 1:256
    
    %Get map with correct pixels "swicthed on"
    bin_idx = 1;
    d = ii-1;
    while d;
        bin_str(bin_idx) = rem(d, 2);
        d = floor(d / 2);
        bin_idx = bin_idx + 1;      
    end
    bin_map(bin_mask) = bin_str;
    
    %Compute normal approximations and nearest pixel
    d_fxi = sum(cx(bin_map));
    d_fyi = sum(cy(bin_map));
    
    d_mag = sqrt(d_fxi^2 + d_fyi^2);
    
    d_fxi = d_fxi / d_mag;
    d_fyi = d_fyi / d_mag;
    
    [~, nearest_px] = min((cx(:) - d_fxi).^2 + (cy(:) - d_fyi).^2);
    
    nearest_x(ii) = cx(nearest_px);
    nearest_y(ii) = cy(nearest_px);
    
    d_fx(ii) = d_fxi;
    d_fy(ii) = d_fyi;
            
    %display(num2str(sum(bin_scores(bin_map(bin_mask)))))
end


