function [rx,ry,rz] = rotlabels(ax,page)
% ROTLABELS Rotate x-, y-, and z-labels to be parallel with the
%
% ROTLABELS rotates the x-, y-, and z-labels of the current 3-Dplot
% to be parallel to the corresponding axis.
%
% ROTLABELS(AX) does the same on the axes with handle ax.
%
% ROTLABELS(AX,PAGE) Any value of PAGE except 0 will force paper
% proportions (as opposed to onscreenproportions).
%
% [RX,RY,RZ] = ROTLABELS(...) will compute the rotation anglesbut
% [R] = ROTLABELS(...) not set them. R = [RX(:) RY(:)RZ(:)].
%
% See also XLABEL, YLABEL, ZLABEL.

% Copyright (c)1995, Erik A. Johnson<johnsone@uxh.cso.uiuc.edu>, 11/21/95

if (nargin<1), ax=gca; end;
if (nargin<2), page=0; elseif (~all(size(page))), page=0; else
page=(page(1)~=0); end;

% could vectorize, but this is simpler for now
if (numel(ax)>1),
    if (nargout==1),
        rx=zeros(numel(ax),3);
    elseif (nargout>1),
        rx=zeros(size(ax));
        ry=zeros(size(ax));
        rz=zeros(size(ax));
    end;
    ax=ax(:);
    for k=1:length(ax);
        if (nargout==1),
            rx(k,:) = rotlabels(ax(k),page);
        elseif (nargout>1),
            [rx(k),ry(k),rz(k)] = rotlabels(ax(k),page);
        else
            rotlabels(ax(k),page);
        end;
    end;
end;

% get axes info
axl = [get(ax,'XLim'); get(ax,'YLim'); get(ax,'ZLim')];
ar = get(ax,'dataAspectRatio');
u = get(ax,'Units');
axposoldunits = get(ax,'Position');
if (page),
    pu = get(get(ax,'Parent'),'PaperUnits');
    set(get(ax,'Parent'),'PaperUnits','points');
    pp = get(get(ax,'Parent'),'PaperPosition');
    set(get(ax,'Parent'),'PaperUnits',pu);
    set(ax,'Units','normalized');
    ap = get(ax,'Position');
    ap = pp.*ap;
else
    set(ax,'Units','pixels');
    ap = get(ax,'Position');
end;
set(ax,'Units',u);
set(ax,'Position',axposoldunits);

% adjust limits for log scale on axes
curxyzlog = [strcmp(get(ax,'XScale'),'log'); ...
             strcmp(get(ax,'YScale'),'log'); ...
             strcmp(get(ax,'ZScale'),'log')];
if (any(curxyzlog)),
    ii = find([curxyzlog;curxyzlog]);
    if (any(axl(ii)<=0)),
        error('ROTLABELS does not support non-positive limits on log-scaled axes.');
    else
        axl(ii) = log10(axl(ii));
    end;
end;

% correct for 'reverse' direction on axes;
curreverse = [strcmp(get(ax,'XDir'),'reverse'); ...
              strcmp(get(ax,'YDir'),'reverse'); ...
              strcmp(get(ax,'ZDir'),'reverse')];
ii = find(curreverse);
if (all(size(ii))),
    axl(ii,[1 2])=-axl(ii,[2 1]);
end;

% correct for aspect ratio
if (~isnan(ar(1))),
    if (ap(3) < ar(1)*ap(4)),
        ap(2) = ap(2) + (ap(4)-ap(3)/ar(1))/2;
        ap(4) = ap(3)/ar(1);
    else
        ap(1) = ap(1) + (ap(3)-ap(4)*ar(1))/2;
        ap(3) = ap(4)*ar(1);
    end;
end;
ap = ap(3:4).';

% correct for 'equal'
% may only want to do this for 2-D views, but seems right for 3-D also
% if (~isnan(ar(2))),
% if ((ap(3)/(axl(1,2)-axl(1,1)))/(ap(4)/(axl(2,2)-axl(2,1)))>ar(2)),
% incr = ap(3)*(axl(2,2)-axl(2,1))/(ap(4)*ar(2)) - (axl(1,2)-axl(1,1));
% axl(1,:) = axl(1,:) + incr/2*[-1 1];
% else,
% incr = ar(2)*(axl(1,2)-axl(1,1))*ap(4)/ap(3) - (axl(2,2)-axl(2,1));
% axl(2,:) = axl(2,:) + incr/2*[-1 1];
% end;
% end;

% transform axes to 2-D space
curT = get(ax,'Xform');
lim = curT*[0 1 0 1 0 1 0 1
            0 0 1 1 0 0 1 1
            0 0 0 0 1 1 1 1
            1 1 1 1 1 1 1 1];
lim = lim(1:2,:)./([1;1]*lim(4,:));
limrange = (max(lim.')-min(lim.')).';

% find right edges
lp = zeros(3);
hc=get(gca,'XLabel'); p=get(hc,'Position'); u=get(hc,'Units');
set(hc,'Units','data'); lp(1,:)=get(hc,'Position');
set(hc,'Units',u); set(hc,'Position',p);
hc=get(gca,'YLabel'); p=get(hc,'Position'); u=get(hc,'Units');
set(hc,'Units','data'); lp(2,:)=get(hc,'Position');
set(hc,'Units',u); set(hc,'Position',p);
hc=get(gca,'ZLabel'); p=get(hc,'Position'); u=get(hc,'Units');
set(hc,'Units','data'); lp(3,:)=get(hc,'Position');
set(hc,'Units',u); set(hc,'Position',p);
ii=find(curxyzlog); if (all(size(ii))), lp(:,ii)=log10(lp(:,ii));
end;
ii=find(curreverse); if (all(size(ii))), lp(:,ii)=-lp(:,ii); end;
[~,ix]=min((lp(1,2)-axl(2,[1 1 2 2])).^2+(lp(1,3)-axl(3,[1 2 1 2])).^2);
[~,iy]=min((lp(2,1)-axl(1,[1 1 2 2])).^2+(lp(2,3)-axl(3,[1 2 1 2])).^2);
[~,iz]=min((lp(3,1)-axl(1,[1 1 2 2])).^2+(lp(3,2)-axl(2,[1 2 1 2])).^2);
ix1 = [1 2; 5 6; 3 4; 7 8];
iy1 = [1 3; 5 7; 2 4; 6 8];
iz1 = [1 5; 3 7; 2 6; 4 8];
ix = [ix1(ix,:); iy1(iy,:); iz1(iz,:)];

% compute 2-D deltas
lim = (lim(:,ix(:,2))-lim(:,ix(:,1)))./limrange(:,[1 1 1]).*ap(:,[1 1 1]);

% compute angles
r = atan2(lim(2,:),lim(1,:))*180/pi;

% adjust for "viewability" -- this may not always be desired.
if (r(1)>90), 
    r(1)=r(1)-180; 
elseif (r(1)<=-90), 
    r(1)=r(1)+180;
end;
if (r(2)>90), 
    r(2)=r(2)-180; 
elseif (r(2)<=-90), r(2)=r(2)+180;
end;
if (r(3)<0),
    r(3)=r(3)+180; 
end;

% do the work
if (nargout==1),
    rx = r;
elseif (nargout>1),
    rx = r(1);
    ry = r(2);
    rz = r(3);
else
    set(get(ax,'XLabel'),'Rotation',r(1));
    set(get(ax,'YLabel'),'Rotation',r(2));
    set(get(ax,'ZLabel'),'Rotation',r(3));
end;

%% code end%%