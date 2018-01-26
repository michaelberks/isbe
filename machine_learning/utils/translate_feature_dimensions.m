function [translated_dims] = translate_feature_dimensions(n1, n2, n3, n4)
%translate_feature_dimensions auxillary tool to translate dimensions
%between different feature representations
%   [out_dims] = translate_feature_dimensions(n1, n2, n3, n4)
%
% Inputs:
%      n1, n2, ... number of features combined factorially
%
%
% Outputs: translated_dims - the translated output dimensions
%
% Example:
%
% Notes: 
%
% See also:
%
% Created: 26-Aug-2015
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
if ~exist('n3', 'var')
    n3 = 1;
end
if ~exist('n3', 'var')
    n4 = 1;
end

n_dims = n1*n2*n3*n4;

translated_dims = zeros(n_dims,1);
out_dim = 1;

for i1 = 0:n1-1
    for i2 = 0:n2-1
        for i3 = 0:n3-1
            for i4 = 0:n4-1
                
                translated_dims(out_dim) = ...
                    i4*n3*n2*n1 + ...
                        i3*n2*n1 + ...
                            i2*n1 + ...
                                i1;
                out_dim = out_dim + 1;
            end
        end
    end
end




