%Circular distributions script

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UNIFORM
% Create a random uniform distribution
sample_u = rand(100000, 1)*2*pi;
figure; rose(sample_u, 50);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VON MISES
% Plot zero-mean Von Mises pdfs with different K values
theta = linspace(0, 2*pi, 200);

p_theta_0pt5 = von_mises_pdf(theta, 0 , 0.5);
p_theta_1 = von_mises_pdf(theta, 0 , 1);
p_theta_2 = von_mises_pdf(theta, 0 , 2);
p_theta_4 = von_mises_pdf(theta, 0 , 4);

figure; hold on;
plot(theta, p_theta_0pt5, 'r');
plot(theta, p_theta_1, 'g');
plot(theta, p_theta_2, 'b');
plot(theta, p_theta_4, 'k');
%%
sample_vm = von_mises_sample(pi, 1, 100000);
figure; hist(sample_vm, theta);
hold on;
plot(theta, 1000*pi*von_mises_pdf(theta, pi, 1), 'r');

%%
%Create some Von-Mises samples
sample_vm4_pi = von_mises_sample(pi, 4, 100000);
figure; rose(sample_vm4_pi, 50);

sample_vm4_0 = von_mises_sample(0, 4, 100000);
figure; rose(sample_vm4_0, 50);

sample_vm0pt5_0 = von_mises_sample(0, .5, 100000);
figure; rose(sample_vm0pt5_0, 50);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Try estimating mu and kappa from a sample drawn from a Von Mises pdf
mu = pi/2; kappa = 4.0;

[sample_vm, cdf] = von_mises_sample(mu, kappa, 1e4);
[mu_hat, kappa_hat] = von_mises_mle(sample_vm);

display(abs(mu - mu_hat));
display(abs(kappa - kappa_hat));
%
theta = linspace(0, 2*pi, 200);
figure; hist(sample_vm, theta);
hold on;
plot(theta, 100*pi*von_mises_pdf(theta, mu_hat, kappa_hat), 'r', 'LineWidth', 2.0);
plot(theta, 100*pi*von_mises_pdf(theta, mu, kappa), 'g:', 'LineWidth', 2.0);
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% CARIOID
theta = linspace(0, 2*pi, 200);

p_theta_1 = cardioid_pdf(theta, pi , .1);
p_theta_2 = cardioid_pdf(theta, pi , .2);
p_theta_4 = cardioid_pdf(theta, pi , .4);

figure; hold on;
plot(theta, p_theta_1, 'g');
plot(theta, p_theta_2, 'b');
plot(theta, p_theta_4, 'k');

%%
sample_car = cardioid_sample(pi, 0.3, 100000);
figure; hist(sample_car, theta);
hold on;
plot(theta, 1000*pi*cardioid_pdf(theta, pi, 0.3), 'r');

%%
sample_car1_pi = cardioid_sample(pi, .1, 100000);
figure; rose(sample_car1_pi, 50);

sample_car1_0 = cardioid_sample(0, .1, 100000);
figure; rose(sample_car1_0, 50);

sample_car4_0 = cardioid_sample(0, .4, 100000);
figure; rose(sample_car4_0, 50);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRAPPED NORMAL
%
theta = linspace(0, 2*pi, 200);

p_theta_1 = wrapped_normal_pdf(theta, pi, .1);
p_theta_2 = wrapped_normal_pdf(theta, pi, .2);
p_theta_4 = wrapped_normal_pdf(theta, pi, .4);

figure; hold on;
plot(theta, p_theta_1, 'g');
plot(theta, p_theta_2, 'b');
plot(theta, p_theta_4, 'k');

%%
[sample_wn] = wrapped_normal_sample(pi, 0.1, 100000);
figure; hist(sample_wn, theta);
hold on;
plot(theta, 1000*pi*wrapped_normal_pdf(theta, pi, 0.1), 'r');

%%
sample_wn1_pi = wrapped_normal_sample(pi, .1, 100000);
figure; rose(sample_wn1_pi, 50);

sample_wn1_0 = wrapped_normal_sample(0, .1, 100000);
figure; rose(sample_wn1_0, 50);

sample_wn4_0 = wrapped_normal_sample(0, .4, 100000);
figure; rose(sample_wn4_0, 50);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRAPPED CAUCHY
theta = linspace(0, 2*pi, 200);

p_theta_1 = wrapped_cauchy_pdf(theta, pi, .1);
p_theta_2 = wrapped_cauchy_pdf(theta, pi, .2);
p_theta_4 = wrapped_cauchy_pdf(theta, pi, .4);

figure; hold on;
plot(theta, p_theta_1, 'g');
plot(theta, p_theta_2, 'b');
plot(theta, p_theta_4, 'k');

%%
[sample_wc] = wrapped_cauchy_sample(pi, 0.1, 100000);
figure; hist(sample_wc, theta);
hold on;
plot(theta, 1000*pi*wrapped_cauchy_pdf(theta, pi, 0.1), 'r');


%%
sample_wc1_pi = wrapped_cauchy_sample(pi, .1, 100000);
figure; rose(sample_wn1_pi, 50);

sample_wc1_0 = wrapped_cauchy_sample(0, .1, 100000);
figure; rose(sample_wn1_0, 50);

sample_wc4_0 = wrapped_cauchy_sample(0, .4, 100000);
figure; rose(sample_wn4_0, 50);


