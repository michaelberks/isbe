% ICASP99A.M
% Plot how step responses vary at level 1 to 4 with 16 shifted steps
% for CWT and DWT.  Animated version.
% Nick Kingsbury, Cambridge University, March 99.
% Modified Dec 2001 to use waitforbuttonpress and Q-shift CWT.

close all
clear all

% Define 16 input steps at staggered positions in columns of x.
n=120;
x=(cumsum([zeros(n,16);diag(ones(16,1));zeros(n,16)]) - 0.5)*1.5;
sc = 1.2;

% DWT of x.
y=dwavxfm(x,4);

% Inverse DWT of each level of y in turn.
yy=y*0;
t=1:(size(y,1)/16);
yy(t,:)=y(t,:);
rdwt = rms(yy);
z0000=dwavixfm(yy,4);
yy=y*0;
t=[1:(size(y,1)/16)] + size(y,1)/16;
yy(t,:)=y(t,:);
rdwt = [rdwt; rms(yy)];
z0001=dwavixfm(yy,4);
yy=y*0;
t=[1:(size(y,1)/8)] + size(y,1)/8;
yy(t,:)=y(t,:);
rdwt = [rdwt; rms(yy)];
z001=dwavixfm(yy,4);
yy=y*0;
t=[1:(size(y,1)/4)] + size(y,1)/4;
yy(t,:)=y(t,:);
rdwt = [rdwt; rms(yy)];
z01=dwavixfm(yy,4);
yy=y*0;
t=[1:(size(y,1)/2)] + size(y,1)/2;
yy(t,:)=y(t,:);
rdwt = [rdwt; rms(yy)];
z1=dwavixfm(yy,4);

zdwt = [x+0.75 z1-1 z01-2 z001-3 z0001-4 z0000-5.75];

% Q-shift DT CWT of x
[C,L] = dtwavedec(x,4,'near_sym_b','qshift_d');

D = C*0; [cm,V] = icwtband(cwtband(C,L,4,'l'),L,4,'l'); D(V) = cm;
z0000 = dtwaverec(D,L,'near_sym_b','qshift_d');
D = C*0; [cm,V] = icwtband(cwtband(C,L,4,'h'),L,4,'h'); D(V) = cm;
z0001 = dtwaverec(D,L,'near_sym_b','qshift_d');
D = C*0; [cm,V] = icwtband(cwtband(C,L,3,'h'),L,3,'h'); D(V) = cm;
z001 = dtwaverec(D,L,'near_sym_b','qshift_d');
D = C*0; [cm,V] = icwtband(cwtband(C,L,2,'h'),L,2,'h'); D(V) = cm;
z01 = dtwaverec(D,L,'near_sym_b','qshift_d');
D = C*0; [cm,V] = icwtband(cwtband(C,L,1,'h'),L,1,'h'); D(V) = cm;
z1 = dtwaverec(D,L,'near_sym_b','qshift_d');

zcwt = [x+0.75 z1-1 z01-2 z001-3 z0001-4 z0000-5.75];

% Select window of output samples for display.
tmax = 40;
tt = -tmax:tmax;

% Draw 1st DWT plot and labels.
setfig(1)
settitle('Shift Dependence')
set(gcf,'DefaultTextFontSize',14,'Color',[1 1 1]);

subplot('position',[0.15 0.1 0.4 0.8]);
plot([-tmax;tmax],[1;1]*[0 -1 -2 -3 -4 -6.5],':k');
axis([-tmax tmax -7 1.6])
axis off
text(0,-7,'Real DWT','horiz','c','FontSize',14)
text(tmax,2,'Reconstructions of a step from coefs at levels 1 to 4 in turn','horiz','c','FontSize',14)
xpos = -42;
text(xpos,0.5,'Input','horiz','r','vert','m');
text(xpos,-0.5,'Wavelets','horiz','r','vert','m');
text(xpos,-1,'Level 1','horiz','r','vert','m');
text(xpos,-2,'Level 2','horiz','r','vert','m');
text(xpos,-3,'Level 3','horiz','r','vert','m');
text(xpos,-4,'Level 4','horiz','r','vert','m');
text(xpos,-5.5,'Scaling fn','horiz','r','vert','m');
text(xpos,-6,'Level 4','horiz','r','vert','m');


st = 1:16:(6*16);

h1(:,1) = line(tt,zdwt(n+8+tt,st),'color','b','LineWidth',2,'erase','none');

pp = -1;

% Draw subsequent DWT plots using h1 line handle.
for i=2:16,
  pause(pp)
  set(h1(:,i-1),'color','c');
  h1(:,i) = line(tt,zdwt(n+8+tt,st+i-1),'color','b','LineWidth',2,'erase','none');
  figure(gcf),drawnow
end

% Do movie of shifting steps with DWT only.
dtime = 0.1;
key = 1;
while key,
   figure(1)
   for i = 16:-1:2,
		ck1 = clock;
      set(h1(:,i),'color','c');
      set(h1(:,i-1),'color','b');
      figure(gcf),drawnow
		while etime(clock,ck1) < dtime, dummy = 1; end
   end
   for i = 1:15,
		ck1 = clock;
      set(h1(:,i),'color','c');
      set(h1(:,i+1),'color','b');
      figure(gcf),drawnow
		while etime(clock,ck1) < dtime, dummy = 1; end
   end
   % [xc,yc] = ginput(1);
   key = waitforbuttonpress; % Press space for repeat, click to move on.
end


% setfig(1)
% set(gcf,'DefaultTextFontSize',12,'Color',[1 1 1]);

% Now add 1st DT CWT plot on right of DWT ones.
subplot('position',[0.58 0.1 0.4 0.8]);
plot([-tmax;tmax],[1;1]*[0 -1 -2 -3 -4 -6.5],':k');
axis([-tmax tmax -7 1.6])
axis off
text(0,-7,'Dual Tree Complex WT','horiz','c','FontSize',14)

st = 1:16:(6*16);

h2(:,1) = line(tt,zcwt(n+8+tt,st),'color','b','LineWidth',2,'erase','none');

% Draw subsequent DT CWT plots using h2 line handle.
for i=2:16,
  pause(pp)
  set(h2(:,i-1),'color','c');
  h2(:,i) = line(tt,zcwt(n+8+tt,st+i-1),'color','b','LineWidth',2,'erase','none');
  figure(gcf),drawnow
end

% Do movie of shifting steps with DWT using h1 and DT CWT using h2.
dtime = 0.1;
key = 1;
while key,
   figure(gcf)
   for i = 16:-1:2,
		ck1 = clock;
      set(h1(:,i),'color','c');
      set(h1(:,i-1),'color','b');
      set(h2(:,i),'color','c');
      set(h2(:,i-1),'color','b');
      figure(gcf),drawnow
		while etime(clock,ck1) < dtime, dummy = 1; end
   end
   for i = 1:15,
		ck1 = clock;
      set(h1(:,i),'color','c');
      set(h1(:,i+1),'color','b');
      set(h2(:,i),'color','c');
      set(h2(:,i+1),'color','b');
      figure(gcf),drawnow
		while etime(clock,ck1) < dtime, dummy = 1; end
   end
   key = waitforbuttonpress; % Press space for repeat, click to move on.
end

return

