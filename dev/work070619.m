%work 19/06/2007

%%
mm = 1;


figure; hist(mass_model.B_shape(mm,:), bins);
cc = cdf('Normal',-1550:100:1550,0, std(mass_model.B_shape(mm,:)));
pp = diff(cc)*185;
hold on;
plot(-1500:100:1500, pp, 'r');
%%
for ii = 1:85
    load(['K:\isbe\dev\annotations\', files(ii).name]);
    mass
end
