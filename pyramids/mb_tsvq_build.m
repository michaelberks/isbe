function [TreeStructure, TreeIndices] = mb_tsvq_build(varargin)

%MB_TSVQ Perform Tree-Structured Vector Quantisation on set of data
%
%
% MB_TSVQ uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%	'Data'
%			- the data to cluster where a row is a vector (i.e. a data point), and
%			the columns are the dimensions of the data.
%
% Optional Arguments:
%	'SimilarityFunction'
%			- 
%	'TotalSize'
%			- a struct of any arguments required by the similarity function
%			being used.
%
%	'SmallClusterSize'
%			- a scalar that describes the minimum size for a cluster; small clusters are
%			discarded. The default value for this is 30 data points.
%	'StoppingPercentage'
%			- 
%
% MB_K_MEANS_CLUSTERING returns:
%	'TreeStructure'
%			- M x (d + 3) array where M is the total number of nodes and d
%			is the dimension of each codeword. Each row corresponds to a
%			node and contains 1) index of 1st child node 2) index of second
%			child node 3) the total error of the codeword (i.e. a measure 
%           of the spread of datapoints around the codeword 4) d-length vector codeword 
%
% References:
%
% pack the args
args = u_packargs(varargin,...
				'0',...
				{'Data'}, ... % the mandatory arguments
				'TotalSize', [],...
                ... the following parameters to be used by Matlab's kmeans
                'Display', 'Notify', ... 
				'SimilarityFunction', 'sqEuclidean', ...
				'Start', 'sample', ...
				'MaxIter', 1000,...
                'Replicates', 1,...
                'EmptyAction', 'singleton',...
                'MaxTreeDepth', 50);
clear varargin;


DATA = args.Data; args = rmfield(args, 'Data');
[n d] = size(DATA);

TreeStructure = zeros(2*(n-1), d+3);
TreeIndices(2*(n-1)).Indices = [];
CurrNode = 0;

args.Indices = 1:n;
args.Node = 0;
args.Depth = 0;
tsvq_recursion(args);

TreeStructure(CurrNode+1:end, :) = [];
TreeIndices(CurrNode+1:end) = [];

display(['Expected nodes = ', num2str(2*(n-1)), ', actual nodes = ', num2str(CurrNode)]);

function tsvq_recursion(recursion_args)

    %First, store node centroid and indices in data structure
    if recursion_args.Node %Don't need to this for root node
        TreeStructure(recursion_args.Node, 3) = 0;
        TreeStructure(recursion_args.Node, 4:end) = recursion_args.Centroid;
        TreeIndices(recursion_args.Node).Indices = recursion_args.Indices;
        
        % Save pointers to child nodes, although these may be overwritten if
        % we turn out to be leaf (this save checking we're not the zero
        % node twice
        TreeStructure(recursion_args.Node, 1) = CurrNode + 1;
        TreeStructure(recursion_args.Node, 2) = CurrNode + 2;
    end

    %Now work out whether we're a leaf node or need to continue recursing
    % 1 of 3 conditions can make us a leaf:
    %  1) We've reached max tree depth
    %  2) We've only got one index left (a true leaf)
    %  3) Clustering returns an empty cluster: i.e. we can't sensibly split
    %  to another level
    
    LeafNode = true;
    
    if recursion_args.Depth >= recursion_args.MaxTreeDepth
        %We're a leaf, save as code -1 in child nodes
        TreeStructure(recursion_args.Node, 1) = -1;
        TreeStructure(recursion_args.Node, 2) = -1;
        
        display(['Max Tree Depth of ', num2str(args.MaxTreeDepth), ' reached' ]);
        
    elseif length(recursion_args.Indices) <= 1
        % We're a 'true' leaf, save as code -2 in child nodes
        
        TreeStructure(recursion_args.Node, 1) = -2;
        TreeStructure(recursion_args.Node, 2) = -2;
        
    elseif length(recursion_args.Indices) > 2
        
        %We've got enough indices to do clustering
        [idx, cluster_means] = kmeans(DATA(recursion_args.Indices, :), 2,...
            'Display', recursion_args.Display,...
            'Distance', recursion_args.SimilarityFunction, 'Start', recursion_args.Start, ...
            'MaxIter', recursion_args.MaxIter, 'Replicates', recursion_args.Replicates, ...
            'EmptyAction', recursion_args.EmptyAction,...
            'Display', 'off');
        
        if all(idx == 1) || all(idx == 2)
            %We're a leaf because we can't split, save as code -3
            TreeStructure(recursion_args.Node, 1) = -3;
            TreeStructure(recursion_args.Node, 2) = -3;
            
        else
            LeafNode = false;
            
            % We can carry on to child nodes
            args_1 = recursion_args;
            args_2 = recursion_args;
            
            args_1.Centroid = cluster_means(1,:);
            args_2.Centroid = cluster_means(2,:);
            
            args_1.Indices = recursion_args.Indices(idx == 1);
            args_2.Indices = recursion_args.Indices(idx == 2);
            
        end
        
    else %length(recursion_args.Indices) == 2, Depth < MaxDepth
        %Don't need to cluster, just split up child indices
        LeafNode = false;
        
        args_1 = recursion_args;
        args_2 = recursion_args;
        
        args_1.Centroid = DATA(recursion_args.Indices(1), :);
        args_2.Centroid = DATA(recursion_args.Indices(2), :);
        
        args_1.Indices = recursion_args.Indices(1);
        args_2.Indices = recursion_args.Indices(2);
        
    end
    
    if ~LeafNode
        %We're not a leaf so...
        
        %Update current node counter
        args_1.Node = CurrNode + 1;
        args_2.Node = CurrNode + 2;
        CurrNode = CurrNode + 2;
        
        %Update tree depth
        args_1.Depth = recursion_args.Depth + 1;
        args_2.Depth = recursion_args.Depth + 1;
        
        %Continue recursing
        tsvq_recursion(args_1);
        tsvq_recursion(args_2);
        
    end
    
%End of recursive function    
end


end