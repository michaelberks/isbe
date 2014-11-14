function [full, unused] = u_packargs(user, strict, varargin)
% U_PACKARGS Fill in optional function arguments as fields in a struct
%   FULL = U_PACKARGS(USER, STRICT, DEFAULT) FULL is a struct which has
%   all the fields which appear in either DEFAULT or USER.  Where a
%   field, F, appears in USER, then FULL.F has the same value as
%   USER.F.  If USER.F does not exist then FULL.F has the same value
%   as DEFAULT.F.  USER can either be a struct, or a cell array.  If
%   the latter then if it is empty then USER is taken to be an empty
%   structure, otherwise it is used to construct a struct thus:
%   STRUCT(USER{:}).  Similarly DEFAULT can either be a struct or a
%   list of (field, value) pairs.  
%
%   STRICT can be true (1 or 'strict'), or false (0 or 'notstrict').
%   If it is false then the behaviour is as described, if true then a
%   field appears in FULL iff it appears in DEFAULT (but its value is
%   as above).  When STRICT is true a warning is given if a field
%   appears in USER but not in DEFAULT (unless there are two output
%   arguments).  This is for use in class constructors, which require
%   that the same set of fields are present and in the same order in
%   each instance.
%
%   Example
%   full = u_packargs(varargin, 'nostrict', 'fig', 1);
%
%   [FULL, UNUSED] = U_PACKARGS(USER, STRICT, DEFAULT) returns
%   the unused arguments from USER in UNUSED.  This will only be
%   non-empty if STRICT is true.
%
%   [FULL, UNUSED] = U_PACKARGS(USER, STRICT, MANDATORY, DEFAULT) Where
%   MANDATORY is a cell array of field names, will produce an error if
%   any of them are not present in USER.  They are always copied to
%   FULL, even if STRICT is true.  In other words they are like fields
%   in DEFAULT except there are no default values associated with them
%   (obviously, because it is an error for them not to be explicitly
%   passed).
%
%   This is useful for filling in optional arguments to functions.
%   DISP_MIP provides an example of its use
%
%   "If you have a procedure with ten parameters, you probably
%   missed some."
%                -- Alan J. Perlis
%
%   "Any sufficiently complicated C or Fortran program contains an
%   ad-hoc, informally-specified bug-ridden slow implementation of
%   half of Common Lisp."
%             -- Phillip Greenspun's Tenth Law of Computing.
%  
%   See also STRUCT, VARARGIN.

  if iscell(varargin{1})
    mandatory = varargin{1};
    varargin = {varargin{2:end}};
  else
    mandatory = {};
  end
  
  if strcmp(strict, 'strict')
    strict = 1;
  elseif strcmp(strict, 'notstrict')
    strict = 0;
  end
  
  args = wrap_cells(varargin);
  default = struct(args{:});
  
  full = default;
  
  % Sort out the USER argument.
  if isempty(user)
    % The user hasn't provided any input.
    full = default;
    unused = [];
    return;
  elseif isstruct(user)
    user = user;
  elseif iscell(user) & length(user) == 1 & isstruct(user{1})
    % Probably a struct passed to a VARARGIN function so it has become
    % wrapped in a cell.
    user = user{1};
  elseif iscell(user) & length(user) == 1 & iscell(user{1})
    args = wrap_cells(user{1});
    user = struct(args{:});
  elseif iscell(user)
    % Presumably (field, value) pairs passed to a VARARGIN
    % function.
    args = wrap_cells(user);
    user = struct(args{:});
  else
    % Might as well try, STRUCT will give a useful error message in
    % any case.
    user = struct(user);
  end
  
  % Ensure that mandatory arguments always appear in the same order.
  for i=1:length(mandatory),
    full = setfield(full, mandatory{i}, []);
  end
  
  % Update FULL from USER.
  names = fieldnames(user);
  mandatory_passed = 0;
  unused = [];
  for i=1:length(names),
    if ismember(mandatory, names{i})
      full = setfield(full, names{i}, getfield(user, names{i}));
      mandatory_passed = mandatory_passed + 1;
    elseif ~strict | isfield(full, names{i})
      full = setfield(full, names{i}, getfield(user, names{i}));
    else
      unused = setfield(unused, names{i}, getfield(user, ...
						   names{i}));
      if nargout < 2
	warning(['Unexpected argument ''' names{i} ''''])
      end
    end
  end
  
  if (mandatory_passed < length(mandatory))
    error('Mandatory argument missing, check the caller and its caller');
  end
  
  
function is = ismember(set, elm)
% Returns true if the strict elm occurs in the cell array set.
  
  for i=1:length(set),
    if strcmp(set{i}, elm)
      is = 1;
      return;
    end
  end
  
  is = 0;
  
function cells = wrap_cells(cells)
% STRUCT interprets cell arrays as a request for an array of
% structs so we need to wrap them up in another cell array with
% only one element.
  
  for i=1:length(cells),
    if iscell(cells{i}) & ~(length(cells{i}) == 1 & iscell(cells{i}))
      cells{i} = {cells{i}};
    end
  end