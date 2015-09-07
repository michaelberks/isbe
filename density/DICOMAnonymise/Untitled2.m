dirName = uigetdir('C:\', 'Choose the image folder');
if (~ischar(dirName))
    return;
end
inUpdatePath=dirName;
nFolder=2;

AllFileName = ls ( sprintf('%s\\',inUpdatePath) );
for i=1:size(AllFileName,1)
    indF=(AllFileName(i,:)~=' ');
    tempInd=find(indF==1);
    indF(1:max(tempInd))=1;
    B=AllFileName(i,indF);
    if ~strcmp(B,'.') && ~strcmp(B,'..') && isdir([inUpdatePath,'\',B]),
       inPathname{nFolder}=[inUpdatePath,'\',B];
       nFolder=nFolder+1;
     %  tmpIn=[inUpdatePath,'\',B];
     %  tmpOut=[outUpdatePath,'\',B]; 
     %  [inPathname, outPathname, nFolder]=checkAllFolder(inPathname, outPathname, nFolder,tmpIn,tmpOut);
    end
end



for ii=2:(nFolder)
    totImage=0;
  
    inFolder=inPathname{ii};
    imageFilename = ls ( sprintf('%s\\*.*',inFolder) );
   
%     k=1;
%     for i=1:size(imageFilename1,1),
%         tmpImageName = sprintf ('%s\\%s',dirName,imageFilename1(i,:));
%         if imageFilename1(i,1)~='.' && ~isdir(tmpImageName)
%           %  if isDICOM(tmpImageName)
%                 imageFilename{k}=imageFilename1(i,:);
%                 k=k+1;
%           %  end
%         end
%     end
    
    %imnum = k-1;
    for i = 3 : size(imageFilename,1),
        name = imageFilename(i,:);
        imageName = sprintf ('%s\\%s',inFolder,name);
        X = dicomread(imageName);
        metadata = dicominfo(imageName);
        
        figure,imshow(X,[]);title([name(1:6),'_',name(8:11),'_',name(12:27),'_',name(29:end)], 'fontsize',10);pause();
    end
        %    [proc_time, pastImnum, subjectID, infoCell,numImage] = AnonmousDICOM (inFolder, [outPathname{i},'\'], totImage, pastImnum, w, subjectID, databaseName, infoCell,numImage);
   
end
ASSURE_MANC_5050TrainingSet_Cancer-00001-LCC121338-20121116-PRO.dcm



