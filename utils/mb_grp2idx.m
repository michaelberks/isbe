function [g,gn,maxgroup] = mb_grp2idx(s)
% GRP2IDX  Create index vector from a grouping variable.
%   [G,GN]=GRP2IDX(S) creates an index vector G from the grouping variable
%   S.  S can be a categorical variable, numeric vector, a character matrix
%   (each row representing a group name), or a cell array of strings stored
%   as a column vector.  The result G is a vector taking integer values
%   from 1 up to the number of unique entries in S.  GN is a cell array of
%   names, so that GN(G) reproduces S (aside from any differences in type).
%
%   See also GRPSTATS, GSCATTER.

%   Copyright 1993-2006 The MathWorks, Inc. 
%   $Revision: 1.4.4.6 $  $Date: 2006/11/11 22:55:11 $

if ischar(s)
   s = cellstr(s);
end

if ~isvector(s)
    error('stats:grp2idx:BadGroup',...
          'Grouping variable must be a vector or a character array.');
end

if isa(s,'categorical')
    g = double(s);
    gn = getlabels(s);
    gn = gn(:);
    maxgroup = length(gn);
    return
end

s = s(:);

%[gn,i,g] = uniquep(s(end:-1:1));           % b=unique group names
[gn,i,g] = unique(s);


if iscell(gn)        % make sure gn is a cell array of strings
    ii = find(strcmp(gn, ''));
    if (isempty(ii))
       ii = find(strcmp(gn, 'NaN'));
    end

    if (~isempty(ii))
       nangrp = ii(1);        % this group should really be NaN
       gn(nangrp,:) = [];     % remove it from the names array
       g(g==nangrp) = NaN;    % set NaN into the group number array
       g = g - (g > nangrp);  % re-number remaining groups
    end
end
maxgroup = length(gn);

% % -----------------------------------------------------------
% function [b,i,j] = uniquep(s)
% % Same as UNIQUE but orders result:
% %    if iscell(s), preserve original order
% %    otherwise use numeric order
% 
% [b,i,j] = unique(s);     % b=unique group names
% i = length(s) + 1 - i; % make sure this is the first instance
% isort = i;  
% if (~iscell(s))  
%    if (any(isnan(b)))  % remove multiple NaNs; put one at the end
%       nans = isnan(b);
%       b = [b(~nans); NaN];
%       x = find(isnan(s));
%       i = [i(~nans); x(1)];
%       j(isnan(s)) = length(b);
%    end
%    isort = b;          % sort based on numeric values
%    if any(isnan(isort))
%       isort(isnan(isort)) = max(isort) + 1;
%    end
% end
% 
% [is, f] = sort(isort); % sort according to the right criterion
% b = b(f,:);
% 
% [fs, ff] = sort(f);    % rearrange j also
% j = ff(j);
% j = j(end:-1:1);
% if (~iscell(b))        % make sure b is a cell array of strings
%    b = cellstr(strjust(num2str(b), 'left'));
% end