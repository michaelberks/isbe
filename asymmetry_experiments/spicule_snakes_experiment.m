%%
%--------------------------------------------------------------------------
%-- Experimental code using snakes to update annotated spicule position
%--------------------------------------------------------------------------

%Load in probability images and spicules
mass_idx = 28;
prob_image = u_load(['M:\chen\data\predict_masses512x512\probability_image', zerostr(mass_idx,3), '.mat']);
spicules = u_load(['M:\chen\data\masses512x512_spicules\mass_spicules', zerostr(mass_idx,3), '.mat']);

if ~isempty(spicules)
    %Set snake parameters
    alpha = 0; %Weight to minimise length (set to zero as don't want to shrink spicule)
    beta = 0.1; %Weight to minimise curvature
    max_delta_x = 3; %Num pixels range horizontally
    resol_x = 1; 
    max_delta_y = 3; %Num pixels range vertically
    resol_y = 1; 
    feat_img = 1 - prob_image; %Feature image = line probability

    %
    for ii = 1:length(spicules)

        %Select i-th spicule
        spicule = spicules{ii};

        %Initialise energies
        e_old = Inf;
        e_new = 1e6;
        thresh = 1e-6;
        iter = 1;

        %Iterate until energy not sufficient reduced from previous iteration
        while e_new+thresh < e_old;

            %Display energy image (i.e line prob) and original spicule in green
            figure; imagesc(feat_img); axis image; colormap(gray(256)); hold on;
            plot(spicule(:,1), spicule(:,2), 'g');
            
            %>>>>>code

            %Store old energy
            e_old = e_new;

            %Apply snake to spicule and compute new energy
            [spicule, e_new] = mb_snake(spicule, alpha, beta, max_delta_x, resol_x, max_delta_y, resol_y, feat_img);

            %Display new spciule in blue/red
            plot(spicule(:,1), spicule(:,2), 'bx');
            plot(spicule(:,1), spicule(:,2), 'r');

            %Print out the new energy and increment the iteration count
            display(['Iteration ', num2str(iter), ': energy = ', num2str(e_new)]);
            iter = iter + 1;
        end
    end
end
%%
% Pre-allocate containers for the normal profiles and the associated
% sampling points
% >>> code
len = size(spicule,1);
norm_width = ?;
normal_p = zeros(len, norm_width);
normal_x = zeros(len, norm_width);
normal_y = zeros(len, norm_width);

%Compute the normal vectors at each point on the inner border
[fx, fy] = gradient(spicule);

%normalise fy
fy = fy ./ [sqrt(sum(fy.^2, 2)), sqrt(sum(fy.^2, 2))];

%Compute normal profiles of the image at every point
for ii = 1:length(fy(:,1)) %= number of rows in skin_air
    
    n1_x = spicule(ii, 1) - norm_width*fy(ii, 2);
    n1_y = spicule(ii, 2) + norm_width*fy(ii, 1);
    n2_x = spicule(ii, 1) + norm_width*fy(ii, 2);
    n2_y = spicule(ii, 2) - norm_width*fy(ii, 1);

    [cx, cy, cp] = improfile(feat_img, [n1_x, n2_x], [n1_y, n2_y], norm_width);
    normal_p(ii, :) = cp';
    normal_x(ii, :) = cx';
    normal_y(ii, :) = cy';

end
normal_p(isnan(normal_p)) = 0;
[spicule,e_new] = mb_snake_normal(spicule, alpha, beta, search_width, resolution, feat_img, normal_x, normal_y);