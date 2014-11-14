function mb_make_new_function(name, use_packargs, args_in, args_out)
%MB_MAKE_NEW_FUNCTION creates a new function with documentation template
%   [] = u_newfunction(function_name,ninput,i1,i2,...,in,o1,o2,...on)
%
%   inputs:
%      name - String giving the name of the new function
%      use_packargs   - Flag if function will use u_packargs
%      args_in - Cell of strings giving the names of the input arguments
%      args_out - Cell of strings giving the names of the output arguments
%
%   outputs:
%   
%   example:
%     to create a function that multiplies x and y to give z:
%     
%     >> mb_make_new_function('mymultiply',0,{'x','y'},{'z'});
%   
%   notes:
%

%Check function name doesn't already exist. If it does, tell user and
%return
if exist(name) %#ok
    display([name, ' already exists.']);
    display('Function not created! Please choose another function name.');
    return;
end

%We can create function...
author = '%% Author: Michael Berks \n';
email = '%% Email : michael.berks@manchester.ac.uk \n';
phone = '%% Phone : +44 (0)161 275 7669 \n';
copyright = '%% Copyright: (C) University of Manchester \n';

if nargin < 2
    use_packargs = 1;
end
if nargin < 3
    args_in = [];
end
if nargin < 4
    args_out = [];
end

filename = [name '.m'];
fid = fopen(filename,'w');

%set list of inputs
if use_packargs
    inputs = ['varargin']; %#ok   
else
    inputs = [];
    for i=1:length(args_in);
       if(i==1)
          inputs = args_in{i};
       else
          inputs = [inputs ', ' args_in{i}]; %#ok
       end
    end
end

%set list of outputs    
outputs = [];
for i=1:length(args_out);
   if(i==1)
      outputs = args_out{i};
   else
      outputs = [outputs ', ' args_out{i}]; %#ok
   end
end

header = ['function [' outputs '] = ' name '(' inputs ')\n'];
summary = ['%%' upper(name) ' *Insert a one line summary here*\n'];
syntax = ['%%   [' outputs '] = ' name '(' inputs ')\n'];

fprintf(fid,header);
fprintf(fid,summary);
fprintf(fid,syntax);
fprintf(fid,'%%\n');

%Print list of input arguments
if use_packargs

    fprintf(fid, ['%% ' upper(name) ' uses the U_PACKARGS interface function\n']);
    fprintf(fid, '%% and so arguments to the function may be passed as name-value pairs\n');
    fprintf(fid, '%% or as a struct with fields with corresponding names. The field names\n');
    fprintf(fid, '%% are defined as:\n');
    fprintf(fid, '%%\n');
    fprintf(fid, '%% Mandatory Arguments:\n');
    fprintf(fid, '%%\n');
    fprintf(fid, '%% Optional Arguments:\n');
    fprintf(fid, '%%\n');
else
    fprintf(fid,'%% Inputs:\n');
    
    for i=1:length(args_in)
       fprintf(fid,['%%      ' args_in{i}]);
       fprintf(fid,' - *Insert description of input variable here*\n');
       fprintf(fid,'%%\n');
    end

    fprintf(fid,'%%\n');
end

%Print list of output arguments
fprintf(fid,'%% Outputs:\n');
for i=1:length(args_out)
   fprintf(fid,['%%      ' args_out{i}]);
   fprintf(fid,' - *Insert description of input variable here*\n');
   fprintf(fid,'%%\n');
end


fprintf(fid,'%%\n');
fprintf(fid,'%% Example:\n');
fprintf(fid,'%%\n');
fprintf(fid,'%% Notes:\n');
fprintf(fid,'%%\n');
fprintf(fid,'%% See also:\n');
fprintf(fid,'%%\n');
fprintf(fid,['%% Created: ' date '\n']);
fprintf(fid,author);
fprintf(fid,email);
fprintf(fid,phone);
fprintf(fid,copyright);

if use_packargs
    fprintf(fid,'\n');
    fprintf(fid,'%% Unpack the arguments:\n');
    fprintf(fid,'args = u_packargs(varargin, 0, {}, {});');
    fprintf(fid,'\n');
    fprintf(fid,'clear varargin;');
end

fclose(fid);

open(name);