function [varargout] = u_load(filename)
%
% U_LOAD - load a MAT file when you don't know the name of the file's variable
%
% Assume you have a MAT file containing data bound to ONE variable, and
% you want to load the data in your code, but don't know how to refer to the
% data (because you don't know the variable name). This function will
% return the data bound to the unknown variable, allowing you to specify
% any variable name for it in your code that you like. An example of when this
% might occur is when two or more pieces of code, which produce the same sort of
% data (usually a struct), use different names for the variable that
% references the data -- when Matlab's SAVE is called, this variable name gets
% saved in the file. If you want to write code that can load data produced
% by the two+ pieces of code, then there is a problem.
%
% An example:
%
% Let 'my-mat-file.mat' be a MAT file containing ONE variable, the name of which
%	is unknown at the time you write the code to load it.
%
% Doing:
%		p = load('my-mat-file.mat')
% would create a structure, p, with one field:
%		p.<unknown_filed_name>
% where <unknown_filed_name> corresponds to the unknown variable name in
% 'my-mat-file.mat'. Since you don't know the name of the variable at the time
% you write the code to load the data, you can't easily write code to cope with
% this situation.
%
% U_LOAD solves this problem:
%		p = u_load('my-mat-file.mat')
%	returns the data, bound to the unknown variable in the file, in the variable
% p. Note that, unlike before, p is not a structure (unless the data in the file
% is actually a structure).
%
% If there is more than one variable in the file, then an error will occur.
% You could probably extend this code to cope with this, but it's beyond my
% needs at the moment.
%
% Extended by Mike Berks to allow multiple variables to be loaded. If there
% are more output arguments than names the function errors. If there are
% more varaibles than output arguments, the function gives a warning and
% loads the first nargout variables into the given names. 
% HOWEVER.... this function is still recommended for .mat files with one 
% variable. Loaded more than one variable requires knowing the names of
% variables, to ensure they are loaded in the correct

% Written by Chris Rose (christopher.rose@stud.man.ac.uk).

% Load the file
S = load(filename);
if isnumeric(S)
    varargout{1} = S;
    return;
end

fnames = fieldnames(S);
if length(fnames) < nargout
	error('The loaded file does not contain enough variables variables!');
elseif length(fnames) > nargout
	warning(['The loaded file contains more variables than output arguments. ',...
        'Only the first ', num2str(nargout)', ' variables loaded.']);
end

for ii = 1:length(fnames)
    % Get the data bound to the unkown (to our caller) variable and return it.
    command = ['S.' fnames{ii}];
    varargout(ii) = {eval(command)}; %#ok
end
