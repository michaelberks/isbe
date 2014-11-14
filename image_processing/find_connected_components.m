function [labels component_sizes] = find_connected_components(bw, connect_size)
%FIND_CONNECTED_COMPONENTS *Insert a one line summary here*
%   [labels, component_sizes] = find_connected_components(bw, connect_size)
%
% Inputs:
%      bw - *Insert description of input variable here*
%
%      connect_size - *Insert description of input variable here*
%
%
% Outputs:
%      labels - *Insert description of input variable here*
%
%      component_sizes - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 28-Oct-2011
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

%copy intial BW image into labels
labels = bw;

%initialise parent_of
parent_of = [];

n_labels = 0;

for ii = 1:size(labels,1)
    for jj = 1:size(labels,2)
        
        if labels(ii,jj)
            
            [n_neighbours neighbours] = get_neighbours(labels, ii, jj, connect_size);
            [min_label max_label] = get_label_range(neighbours);
            
            if min_label == -1
                %pixel has no connections, create a new label
                n_labels = n_labels + 1;
                labels(ii,jj) = n_labels;
                parent_of(end+1) = 0; %#ok
                
            elseif min_label == max_label
                %pixel connected to one known component, assign to this
                %component
                labels(ii,jj) = min_label;
                
            else
                %pixel connects 2 different components
                labels(ii,jj) = min_label;
                
                roots = neighbours;
                [parent_of roots] = get_roots(parent_of, roots, n_neighbours);
                
                %get lowest root index
                min_root = min(roots);
                
                %make this root the parent of all the others
                for nn = 1:n_neighbours
                    if roots(nn) ~= min_root
                        parent_of(roots(ii) - 1) = min_root;
                    end
                end
            end
        end
    end
end

true_labels = zeros(length(parent_of),1);
n_components = 0;
for ii = 1:length(parent_of)
    if ~parent_of(ii)
        %root node
        n_components = n_components + 1;
        true_labels(ii) = n_components;
    
    else
        %track back up the tree to find this node's root
        
        %check neither this node nor its parent is a root node
        while (parent_of(ii) > 0) && (parent_of( parent_of(ii) - 1) > 0)
            parent_of(ii) = parent_of( parent_of(ii) - 1 );
        end
        
        true_labels(ii) = true_labels( parent_of(ii) - 1 );
    end
end

component_sizes = zeros(n_components, 1);

%Update label image and record component sizes
for ii = 1:size(labels,1)
    for jj = 1:size(labels,2)
        
        if labels(ii,jj)
            
            labels(ii,jj) = true_labels( labels(ii,jj) - 1 );
            component_sizes(labels(i,j)) = component_sizes(labels(i,j)) + 1;
        end
    end
end
                
            



function [n_neighbours neighbours] = get_neighbours(labels, ii, jj, connect_size)

n_neighbours = 0;
neighbours = [];
%If not on top edge, add pixel above
if (ii > 1) && (labels(ii-1,jj) > 0)
    n_neighbours = n_neighbours + 1;
    neighbours(n_neighbours) = labels(ii-1,jj);   
end

%If not on left edge, add pixel to left
if (jj > 1) && (labels(ii,jj-1) > 0)
    n_neighbours = n_neighbours + 1;
    neighbours(n_neighbours) = labels(ii,jj-1);   
end

%If using 8-neighbourhood, add pixel above left/right
if (connect_size == 8) && (ii > 1)
    
    if (jj > 1) && (labels(ii-1,jj-1) > 0)
        n_neighbours = n_neighbours + 1;
        neighbours(n_neighbours) = labels(ii-1,jj-1);   
    end

    if (jj < size(labels,2)) && (labels(ii-1,jj+11) > 0)
        n_neighbours = n_neighbours + 1;
        neighbours(n_neighbours) = labels(ii-1,jj+1);   
    end
end

function [parent_of roots] = get_roots(parent_of, roots, n_neighbours)

for ii = 1:n_neighbours
    while parent_of( roots(ii) - 1) > 0
        roots(ii) = parent_of( roots(ii) - 1 );
    end
end

function [min_label max_label] = get_label_range(neighbours)

if isempty(neighbours)
    min_label = -1;
    max_label = -1;
else
    min_label = min(neighbours);
    max_label = max(neighbours);
end
    