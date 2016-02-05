function [ContentList] = ReadContentOfDirectory(DataPath,keyword)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% function ReadContentOfDirectory
% 
% The function reads the content of a directory and sorts it in order. 
% The function is mainly code from the pda_ReadImageDirectory from the
% LC/MS project
%
% syntax: [ContentList] = ReadContentOfDirectory(DataPath,keyword)
% 
% Inputs: 
%        DataPath: The path of the directory
%        keyword: It helps in retrieving content with a specific keyword
%
% Output:
%        ContentList: The content (names of the files) of the directory 
%                   sorted in a list  
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
flag=0;

% Reads the contents of a directory and sorts them in order
% disp('Reading the content of the Directory'); 
if nargin < 2
    List = dir(DataPath);
    
    % The first two will be './' and '../' so we can discard those
    % JG It doesn't seem to be that simple. There are cases where ./ and ../
    %appear at the end of the list.  code modified to make sure they don't
    %finish up in the list.
    ImageCount = 1;
    for i = 1:length(List)
        Name = List(i).name;
        if ~strcmp(Name,'Thumbs.db') && ~strcmp(Name,'.') && ~strcmp(Name,'..');
            ContentList{ImageCount} = Name;
            ImageCount = ImageCount + 1;
            flag=1;
        end
    %Danny's version
%     for i = 3:length(List)
%         Name = List(i).name;
%         %if Name(end-1:end) ~= 'db';
%         if ~strcmp(Name,'Thumbs.db');
%             ContentList{ImageCount} = Name;
%             ImageCount = ImageCount + 1;
%         end
    end
    %if there are no files in the folders return an empty list (EM)
    if flag==0
        ContentList=[];
    end
else
    List = dir(strcat(DataPath,'/*',keyword,'*'));
    ImageCount = 1;
    for i = 1:length(List)
    Name = List(i).name;
    %if Name(end-1:end) ~= 'db';
        if ~strcmp(Name,'Thumbs.db');
            ContentList{ImageCount} = Name;
            ImageCount = ImageCount + 1;
            flag=1;
        end
    end
    %if there are no files in the folders return an empty list (EM)
    if flag==0
        ContentList=[];
    end
end