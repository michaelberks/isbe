function directions = sample_direction(I_prob, x, y, direction0, ...
			  						   N, n_bins, step_length, theta_max, theta_min)

min_args = 4;
if nargin<min_args+1, N = 1; end
if nargin<min_args+2, n_bins = 36; end
if nargin<min_args+3, step_length = 4; end
if nargin<min_args+4, theta_max = pi/2; end
if nargin<min_args+5, theta_min = -theta_max; end

% get current heading
theta0 = angle(direction0);

% compute vector of potential new headings
theta_vec = theta0 - linspace(theta_min, theta_max, n_bins);

% convert these to pixel positions
x_vec = x + step_length*cos(theta_vec);
y_vec = y - step_length*sin(theta_vec);

% sample image at these positions
p_vec = interp2(I_prob, x_vec, y_vec, '*linear');

% figure(2); hold on;
% 	plot(x_vec,y_vec,'y.');
% 	plot(x_vec(1),y_vec(1),'g.');

% sample new directions
inds = sample_from_pvec(p_vec);

% compute new headings and return
directions = theta_vec(inds);
