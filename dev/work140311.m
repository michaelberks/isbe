load('C:\isbe\nailfold\data\rsa_study\master_set\apex_maps\set12g_half_296655\miccai_maxima\results\cooccuurence.mat');

n = 100;
b = 1000;
m = 10000;
km = zeros(100,3);

for ii = 1:n
    [km(ii,1) km(ii,2) km(ii,3)] = miccai_kappa_stats(co_occurrence_t, b, m*(ii-1));
end
figure;
plot(m*(0:n-1)', km);

title({'Kappa agreement for marking vessels as a function of the ';...
    'structural zero component of the co-occurrence matrix (C_{1,1,1})',});
xlabel('C_{1,1,1}');
ylabel('Kappa agreement - \kappa_m(i,j)');
legend({'\kappa_m(1,2)', '\kappa_m(1,3)', '\kappa_m(2,3)'}, 'location', 'southeast');
exportfig('M:\nailfold\weekly_presentations\figures\kappa_agreement_limits.png');
%%
co = co_occurrence_t;
co(1,1,1) = co(1,1,1) + 1e6;
[km_l1 km_l2 km_l3] = miccai_kappa_stats(co, b);

% hold on;
% plot([1e5 1e6], [km(end,1) km_l1], 'b');
% plot([1e5 1e6], [km(end,2) km_l2], 'g');
% plot([1e5 1e6], [km(end,3) km_l3], 'r');
%%
co = sum(sum(co_occurrence_t,4),3);
a = co(1,1) + 1e6;
b = co(1,2) + co(1,3);
c = co(2,1) + co(3,1);
d = co(2,2) + co(2,3) + co(3,2) + co(3,3);

%a = linspace(0, 1e6, 1000);
k = 2*(a*d - b*c) ./ (2*a*d + b*b + c*c + a*(b+c) + d*(b+c));
k_lim = 2*d/(2*d + b + c);
figure; plot(a, k); hold on;
plot([0 1e6], [k_lim k_lim], 'r--');


