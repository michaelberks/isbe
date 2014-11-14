function [Ix, Iy, It] = image_derivatives(varargin)
% This is now a wrapper function for differentiate_stack.

[Ix, Iy, It] = differentiate_stack(varargin{:});
