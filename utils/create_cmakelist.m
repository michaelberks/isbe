function [] = create_cmakelist(src_dir, project_name, varargin)
%CREATE_CMAKELIST *Insert a one line summary here*
%   [] = create_cmakelist(src_dir, project_name)
%
% Inputs:
%      src_dir - directory of source code
%
%      project_name - name of project to be generate in CMake
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
% Created: 05-Apr-2017
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester

args = u_packargs(varargin, '0',...
    'libraries', [],...
    'executables', [],...
    'use_qt4', 0,...
    'use_qt4_wrap_ui', 0,...
    'use_qt4_wrap_cpp', 0,...
    'use_qt4_add_resources', 0,...
    'cmake_min_ver', '2.8',...
    'author', '# Author: Michael Berks \n',...
    'email', '# Email : michael.berks@manchester.ac.uk \n',...
    'phone', '# Phone : +44 (0)161 275 7669 \n',...
    'copyright', '# Copyright: (C) University of Manchester \n',...
    'overwrite', 0);
clear varargin;

filename = fullfile(src_dir, 'CMakeLists.txt');
if exist(filename, 'file') && ~args.overwrite
    display([filename, ' already exists, and overwrite is false, new file CMakeLists file not created.']);
    return;
end

fid = fopen(filename,'w');

%Some header info for file
fprintf(fid,['# Created: ' date '\n']);
fprintf(fid, args.author);
fprintf(fid, args.email);
fprintf(fid, args.phone);
fprintf(fid, args.copyright);

%Write project name
fprintf(fid,'\n');
fprintf(fid,['project(' project_name ')\n']);

%Set required cmake version
if ~isempty(args.cmake_min_ver)
    fprintf(fid,'\n');
    fprintf(fid,['cmake_minimum_required(VERSION ' args.cmake_min_ver ')\n']);
end

if args.use_qt4
    fprintf(fid,'\n');
    fprintf(fid,'# Add support for Qt4\n');
    fprintf(fid,'INCLUDE(${CMAKE_ROOT}/Modules/FindQt4.cmake)\n');
    fprintf(fid,'INCLUDE(${QT_USE_FILE}\n)');
end

if args.use_qt4_wrap_ui
    fprintf(fid,'\n');
    fprintf(fid,'# Convert these .ui files to .h\n');
    fprintf(fid,['QT4_WRAP_UI( ' project_name '_ui_h_files \n']);
    fprintf(fid,')\n');
end	

if args.use_qt4_wrap_cpp
    fprintf(fid,'\n');
    fprintf(fid,'# Convert these .h files to .cxx\n');
    fprintf(fid,['QT4_WRAP_CPP( ' project_name '_moc_files \n']);
    fprintf(fid,')\n');
end	

if args.use_qt4_add_resources
    fprintf(fid,'\n');
    fprintf(fid,'# Convert these resource files\n');
    fprintf(fid,['QT4_ADD_RESOURCES( ' project_name '_res_files \n']);
    fprintf(fid,')\n');
end

header_file_list = dir(fullfile(src_dir, '*.h'));
cxx_file_list = [...
    dir(fullfile(src_dir, '*.cxx')) ;
    dir(fullfile(src_dir, '*.cpp')) ;
    dir(fullfile(src_dir, '*.c')) ];

header_files = {header_file_list(:).name}';
cxx_files = {cxx_file_list(:).name}';

[~, header_filenames, header_file_ext] = ...
    cellfun(@fileparts, header_files, 'UniformOutput', false);

[~, cxx_filenames, cxx_file_ext] = ...
    cellfun(@fileparts, cxx_files, 'UniformOutput', false);

all_filenames = union(header_filenames, cxx_filenames);

fprintf(fid,'\n');
fprintf(fid,['SET(' project_name '_sources\n']);

for i_f = 1:length(all_filenames)
    name_i = all_filenames{i_f};
    
    cidx = find(strcmpi(cxx_filenames, name_i));
    if ~isempty(cidx)
        fprintf(fid,['    ' name_i cxx_file_ext{cidx}]);
    end
    
    hidx = find(strcmpi(header_filenames, name_i));
    if ~isempty(hidx)
        fprintf(fid,['    ' name_i header_file_ext{hidx}]);
    end
    
    fprintf(fid,'\n');
end
fprintf(fid,')\n');


for i_exe = 1:length(args.executables)
    fprintf(fid,['ADD_EXECUTABLE( ' args.executables{i_exe} '\n']);
    
    fprintf(fid,['    ${' args.executables{i_exe} '_sources} \n']);
    if args.use_qt4_wrap_ui
        fprintf(fid,['    ${' args.executables{i_exe} '_ui_h_files} \n']);
    end	
    if args.use_qt4_wrap_cpp
        fprintf(fid,['    ${' args.executables{i_exe} '_moc_files} \n']); 
    end	
    if args.use_qt4_add_resources
        fprintf(fid,['    ${' args.executables{i_exe} '_res_files} \n']);
    end
	
    fprintf(fid,')\n');
    fprintf(fid,['TARGET_LINK_LIBRARIES( ' args.executables{i_exe} ' )\n']);
end

for i_exe = 1:length(args.libraries)
    fprintf(fid,['ADD_LIBRARY( ' args.executables{i_exe} '\n']);
	fprintf(fid,['    ${' args.executables{i_exe} '_sources} \n']);
    if args.use_qt4_wrap_ui
        fprintf(fid,['    ${' args.executables{i_exe} '_ui_h_files} \n']);
    end	
    if args.use_qt4_wrap_cpp
        fprintf(fid,['    ${' args.executables{i_exe} '_moc_files} \n']); 
    end	
    if args.use_qt4_add_resources
        fprintf(fid,['    ${' args.executables{i_exe} '_res_files} \n']);
    end
    fprintf(fid,')\n');
end

%we're done!
fclose(fid);
open(filename);
