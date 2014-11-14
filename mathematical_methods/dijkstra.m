function [path_list, path_lengths] = dijkstra(first, connections)
%DIJKSTRA Dijkstra's algorithm for shortest paths through a graph
%  [PATHS, PATH_LENGTHS] = DIJKSTRA(FIRST, CONNECTIONS)
%  
%  inputs:
%           FIRST - Scalar, first node
%     CONNECTIONS - (NxN) each element the length between nodes
%                   (+inf if the nodes are not directly connected).
%  
%  outputs:
%    PATHS - (1xN) cell array of paths such that PATH{i} is a vector
%           of nodes giving the shortest path from FIRST to i.
%    PATH_LENGTHS - (1xN) vector of the lengths of each path.
% 
%  N is the number of junctions.  Dijkstra's algorithm is O(N^2) in
%  principle but this implementation is not.
%
%  Example:
%    A graph with four nodes with d(1,2)=8, d(1,3)=5, d(2,3)=1, d(1,4)=7
%    paths = dijkstra(1, [0 8 5 7; 8 0 1 2; 5 1 0 inf; 7 2 inf 0])
%    
% 
%  References:
%    Erwin Kreyszig, Advanced Engineering Mathematics, 7th ed., p.
%    1118.
  
  % Number of junctions (vertices)
  nj = size(connections, 1);
  
  % ls(i, j) is the length of the edge between vertices i and j.
  ls = connections;
  
  % The first vertex
  v1 = first;

  % tl(i)=i, iff the path to i has not yet been found
  % tl(i)=0, otherwise.
  % In other words, find(tl) gives a list of vertices yet to be
  % processed.
  tl = 1:nj;
  tl(v1) = 0;
  
  % For finite L_tmp(i), L_tmp(i) is the length of the shortest path
  % yet found between v1 and i.  If L_tmp(i) is infinite then we have
  % found the overall shortest path to vertex i, the length of which
  % is given by L(i).
  L_tmp = ls(v1, :);
  L_tmp(v1) = inf;
  
  % L(i) is the length of the shortest path between v1 and i, or inf
  % if the path has not yet been found.
  L = repmat(inf, 1, nj);
  L(v1) = 0;
  
  % For all i,
  % tl(i)=0 <=> L_tmp(i)=inf <=> L(i)<inf
  % tl(i)=i <=> L_tmp(i)<inf <=> L(i)=inf
  
  % path_list{i} will be a list of vertices defining the shortest path
  % between v1 and i, or {} if there is no such path.
  path_list = cell(nj, 1);
  
  while (1),
    % If there are no more possible paths, but we have not yet visited
    % all the vertices then there must be some disconnected
    % sub-graphs.  We will return an empty path for them below.
    if ~any(isfinite(L_tmp))
      break;
    end

    % Find the shortest remaining path...
    [mn, k] = min(L_tmp);
    k = k(1);
    % and save it as the best path between v1 and k.
    L(k) = L_tmp(k);
    % Remove k from the list of unprocessed vertices.
    L_tmp(k) = inf;
    tl(k) = 0;
    
    % If we have visited all vertices then we are done.
    if ~any(tl)
      break
    end
    
    % Update all the lengths for the remaining vertices.
    for j=find(tl),
      % If the path through vertex k to vertex j is better than the
      % current best path...
      if L(k) + ls(j, k) < L_tmp(j)
	% ...replace the current path
	L_tmp(j) = L(k) + ls(j, k);
	path_list{j} = [path_list{k}, k];
      end
    end
  end
  
  % Add the final node to each path
  for j=1:length(path_list),
    path_list{j} = [first, path_list{j}, j];
  end
  % Must set this because otherwise path_list{first} would be set to
  % [first first] by the loop above.
  path_list{first} = first;

  % There is no path from v1 to the vertices still in tl.
  for j=find(tl),
    path_list{j} = {};
  end
  
  path_lengths = L;

