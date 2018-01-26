function [patch_handles] = voxel(v_xyz, dxyz, face_color, alpha, edge_color)

%VOXEL function to draw a 3-D voxel in a 3-D plot
%
%Usage
%   voxel(start,size,color,alpha);
%
%   will draw a voxel at 'start' of size 'size' of color 'color' and
%   transparency alpha (1 for opaque, 0 for transparent)
%   Default size is 1
%   Default color is blue
%   Default alpha value is 1
%
%   start is a three element vector [x,y,z]
%   size the a three element vector [dx,dy,dz]
%   color is a character string to specify color 
%       (type 'help plot' to see list of valid colors)
%   
%
%   voxel([2 3 4],[1 2 3],'r',0.7);
%   axis([0 10 0 10 0 10]);
%

%   Suresh Joel Apr 15,2003
%           Updated Feb 25, 2004

%Set defaults, then swap out as args are supplied   
if nargin < 5
    edge_color = 'none';
end
if nargin < 4
    alpha=1;
end
if nargin < 3
    face_color='b';
end
if nargin < 2
    dxyz=eye(3);
end
if nargin < 1
    disp('Too few arguments for voxel');
    return;
end

x = [...
    0 0 0;
    0 0 1;
    0 1 0;
    0 1 1;
    1 0 0;
    1 0 1;
    1 1 0;
    1 1 1] * dxyz + repmat(v_xyz, 8, 1);

patch_handles = zeros(1,6);
curr_face = 1;

for n=1:3,
    if n==3,
        x=sortrows(x,[n,1]);
    else
        x=sortrows(x,[n n+1]);
    end;
    temp=x(3,:);
    x(3,:)=x(4,:);
    x(4,:)=temp;
    h=patch(x(1:4,1),x(1:4,2),x(1:4,3),face_color, 'EdgeColor', edge_color);
    patch_handles(curr_face) = h;
    curr_face = curr_face+1;
    
    set(h,'FaceAlpha',alpha);
    temp=x(7,:);
    x(7,:)=x(8,:);
    x(8,:)=temp;
    h=patch(x(5:8,1),x(5:8,2),x(5:8,3),face_color, 'EdgeColor', edge_color);
    set(h,'FaceAlpha',alpha);
    patch_handles(curr_face) = h;
    curr_face = curr_face+1;
end;