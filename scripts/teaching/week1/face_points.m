load('P:\MATLAB\data\face_points.mat');

%face_x and face_y are 68x1 vectors, specifying the points of a face

%Create a new point
pt_xy = [70 60];

%Which point is [pt_xy(:,1) pt_xy(:,2)] closest to? 
dists = sqrt((face_x-pt_xy(:,1)).^2 + (face_y-pt_xy(:,2)).^2); %Using elementwise square

%This is an example of using 2 outputs from a function - we'll return to
%this in the next lesson
[min_dist, min_i] = min(dists);

x_nearest = face_x(min_i);
y_nearest = face_y(min_i);
%--------------------------------------------------

%Now lets rotate the points by theta degrees
theta = 30;

%Create the rotation matrix - cosd, sind apply cosine and sine to degree
%inputs (use cos, sin for radians)
R = [cosd(theta) sind(theta); -sind(theta) cosd(theta)];
xy_rot = [face_x face_y]*R; %Using matrix multiplication, xy_rot is 68x2
x_rot = xy_rot(:,1);
y_rot = xy_rot(:,2);

%--------------------------------------------------
%Now lets visually display our output - we'll look at plotting commands in
%more detail in lesson 7

figure; %Create a new window;
axis equal; %Shortcut command to make the face_x and face_y axes equal in the plot;
hold on; %Shortcut command so that subsequent plots add to the axes rather than overwriting them
plot(face_x, face_y, 'rx'); %Show each [face_x,face_y] point with a red cross
plot(face_x, face_y, 'b-'); %Show the points joined by blue lines (like a dot-the-dot puzzle!)
plot(pt_xy(:,1), pt_xy(:,2), 'gx'); %Show [pt_xy(:,1) pt_xy(:,2)] with a green cross
plot(x_nearest, y_nearest, 'go'); %draw a green circle round the point closest to [pt_xy(:,1) pt_xy(:,2)]
plot([pt_xy(:,1) x_nearest], [pt_xy(:,2) y_nearest], 'g-') %Connect the above with a green line

%Create a new plot and draw the rotated face
figure;
axis equal;
hold on;
plot(x_rot, y_rot); %Note that 'b-' is the default setting

%PS. If you're wondering what all this has to do with image processing,
%these points are used to buildmodels capable of face-tracking - as
%demonstrated in this video by my colleague Phil:
% http://www.youtube.com/watch?v=5TDO9ok4sWI
