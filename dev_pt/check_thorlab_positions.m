clc;

dat = dlmread('u:\tmp\positions.txt');
t = (dat(:,1) - dat(1,1)) / 1000.0;
p = dat(:,2);

figure(1); clf; 
subplot(3,1,1);
    plot(t, p);
    axis([0,max(t),0,10.0]);

dt = median(t(3:end) - t(1:end-2));

dp = p(3:end) - p(1:end-2);
v = dp / dt;
maxv = 1.2;
subplot(3,1,2); hold on;
    plot(t(2:end-1,1), v, 'b-');
    plot(t(2:end-1,1), dat(2:end-1,3), 'r-');
    axis([0,max(t),-maxv,maxv]);

dv = v(3:end) - v(1:end-2);
a = dv / dt;
subplot(3,1,3); hold on;
    plot(t(3:end-2,1), a, 'b-');
    plot(t(3:end-2,1), dat(3:end-2,4), 'r-');
    axis([0,max(t),-2.0,2.0]);
