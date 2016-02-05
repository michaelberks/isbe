function ImageNamePathMat=Record_Image_Names_and_Paths(search_arg, verbose)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Function Record_Image_Names_and_Paths
% 
% The function records the images and paths of
% the images in the subfolders of a folder.
%  
% syntax: ImageNamesPathMat=Record_Image_Names_and_Paths(search_arg,verbose)
% 
% intput: 
%        verbose: if 1, the function reports the findings in the workspace
%        search_arg: text (e.g. 'RAW.dcm') that specifies the search
%        criterion
% 
% Emmanouil Moschidis 10/11/2013
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%List the subdirectories of a directory in a cell (folderNames)
currentPath=cd();                       % record the current directory
dlist = dir(currentPath);               % list all the directories
subvec = [dlist(:).isdir];              % make sure there are only directories (subvec is a logical vector)
folderNames = {dlist(subvec).name}';    % list all the directory names 
folderNames(ismember(folderNames,{'.','..'})) = []; %the first two directories ('.' and '..') are truncated
counter=0;
%loop through the folders 
for i=1:length(folderNames) 
    cd (folderNames{i,1}); %move to the appropriate subdirectory
    newpath=cd();          %get the full path
    ListOfImages=ReadContentOfDirectory(newpath,search_arg); %get a list of the images according to the search argument
    %loop through the images of the directory
    for cntr=1:length(ListOfImages)
        %info = dicominfo(ListOfImages{1,cntr});                 %get the DICOM header
        if verbose==1; 
            disp(['Acessing image number: ',num2str(counter+1), ' with filename: ' , ListOfImages{1,cntr}]);
        end;  
        counter=counter+1;
        ImageNamePathMat{counter,1}=ListOfImages{1,cntr};
        ImageNamePathMat{counter,2}=cd();
    end
    clear ListOfImages;
    cd(currentPath);       %return to the parent directory 
end


