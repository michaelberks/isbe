function [tout rout] = polar_bar(bin_counts, bins)
%POLAR_BAR *Insert a one line summary here*
%   [] = polar_bar(bin_counts, bins)
%
% Inputs:
%      bin_counts - *Insert description of input variable here*
%
%      bins - *Insert description of input variable here*
%
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 28-Oct-2010
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester

%Check if bin counts is a vector of array
if min(size(bin_counts))==1, % For vectors arrange as a single column vec
    bin_counts = bin_counts(:); 
end

%See if bins is scalar (in which case it specifies the number of bins) or
%is a vector specifying the bin centres themselves
if length(bins)==1
    %If scalar, create bins uniformly spread of the unit circle
    bins = (0:bins-1)*2*pi/bins + pi/bins;
else
    %If a vector, divide modulo 2*pi to ensure all bin centres lie within
    %the unit circle
    bins = sort(rem(bins(:)',2*pi));
end

%Compute the bin edges spaced evenly between each bin (and wrapping round
%at 2*pi)
edges = sort(rem([(bins(2:end)+bins(1:end-1))/2 (bins(end)+bins(1)+2*pi)/2],2*pi));
edges = [edges edges(1)+2*pi];

% Form rho,theta values for histogram triangle given the bin edges and bin
% counts

%Rho values
[m,n] = size(bin_counts);
mm = 4*m;
rho = zeros(mm,n);
rho(2:4:mm,:) = bin_counts;
rho(3:4:mm,:) = bin_counts;

%theta values
theta = zeros(mm,1);
theta(2:4:mm) = edges(1:m);
theta(3:4:mm) = edges(2:m+1);

%If we're not returning rho and theta, plot the histogram in polar form and
%return plot handle
if nargout < 2
    h = polar(theta, rho);
    if nargout==1
        tout = h;
    end
else
    %Return rho and theta
    tout = theta; rout = rho;
end