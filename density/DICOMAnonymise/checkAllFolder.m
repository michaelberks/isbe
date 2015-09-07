function [inPathname, outPathname, nFolder]=checkAllFolder(inPathname, outPathname, nFolder,inUpdatePath,outUpdatePath)

AllFileName = ls ( sprintf('%s\\',inUpdatePath) );
for i=1:size(AllFileName,1)
    indF=(AllFileName(i,:)~=' ');
    tempInd=find(indF==1);
    indF(min(tempInd):max(tempInd))=1;
    B=AllFileName(i,indF);
    if ~strcmp(B,'.') && ~strcmp(B,'..') && isdir([inUpdatePath,'\',B]),
       inPathname{nFolder}=[inUpdatePath,'\',B];
       outPathname{nFolder}=[outUpdatePath,'\',B]; 
      % mkdir(outUpdatePath,B);
       nFolder=nFolder+1;
       tmpIn=[inUpdatePath,'\',B];
       tmpOut=[outUpdatePath,'\',B]; 
       [inPathname, outPathname, nFolder]=checkAllFolder(inPathname, outPathname, nFolder,tmpIn,tmpOut);

    end
end
