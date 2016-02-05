function [ I1 J1 I2 J2 ] = polar2borderpoints(imgsize,line_rho,line_theta)
% POLAR2BORDERPOINTS returns the end points of a line given by rho and
% theta in an image of size imgsize.
%
% Andreas Eilschou, 2010

% The line intersects with two of the four borders.
line = zeros(2,2);
    line_point = 1;
    J = 1;
    I = (line_rho - (J-0.5) * cos(line_theta)) / sin(line_theta) + 0.5;
    if (I >= 1 && I <= imgsize(1))
        line(line_point,:) = [I,J];
        line_point = line_point + 1;
    end
    J = imgsize(2);
    I = (line_rho - (J-0.5) * cos(line_theta)) / sin(line_theta) + 0.5;
    if (I >= 1 && I <= imgsize(1))
        line(line_point,:) = [I,J];
        line_point = line_point + 1;
    end
    I = 1;
    J = (line_rho - (I-0.5) * sin(line_theta)) / cos(line_theta) + 0.5;
    if (line_point <= 2 && J >= 1 && J <= imgsize(2))
        line(line_point,:) = [I,J];
        line_point = line_point + 1;
    end
    I = imgsize(1);
    J = (line_rho - (I-0.5) * sin(line_theta)) / cos(line_theta) + 0.5;
    if (line_point <= 2 && J >= 1 && J <= imgsize(2))
        line(line_point,:) = [I,J];
    end
    I1 = line(1,1);
    J1 = line(1,2);
    I2 = line(2,1);
    J2 = line(2,2);
end