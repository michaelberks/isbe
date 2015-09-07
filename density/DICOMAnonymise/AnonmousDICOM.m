function [e, pastImnum, infoCell,numImage] = AnonmousDICOM (dirName, savePathname, totImage, pastImnum, w, subjectID, databaseName,infoCell,numImage)
warning off;

imageFilename1 = ls ( sprintf('%s\\*.*',dirName) );
if isempty(imageFilename1)
    e = 0;
    return;
end

k=1;
for i=1:size(imageFilename1,1),
    tmpImageName = sprintf ('%s\\%s',dirName,imageFilename1(i,:));
    if imageFilename1(i,1)~='.' && ~isdir(tmpImageName)
        if isDICOM(tmpImageName)
            imageFilename{k}=imageFilename1(i,:);
            k=k+1;
        end
    end
end

if ~exist('imageFilename')
    e = 0;
    return;
end
imnum = k-1;
t = clock;
for i = 1 : imnum
    name = imageFilename{i};
    imageName = sprintf ('%s\\%s',dirName,name);
    X = dicomread(imageName);
    metadata = dicominfo(imageName);
    if exist('metadata', 'var')
        
        iter=5-length(num2str(subjectID));
        subStr=num2str(subjectID);
        for j=1:iter,
            subStr=['0',subStr];
        end
        if isfield(metadata,'StudyDate')
            scanDate=metadata.StudyDate;
        else
            scanDate='--';
        end
        if isfield(metadata,'ViewPosition')
            View=metadata.ViewPosition;
        else
            View='--';
        end
        if isfield(metadata,'ImageLaterality')
            LR=metadata.ImageLaterality;
        else
            LR='--';
        end
        if isfield(metadata,'AcquisitionTime')
            sTime=metadata.AcquisitionTime(1:6);
        else
            sTime='--';
        end
        if isfield(metadata,'PresentationIntentType')
            if strcmp(metadata.PresentationIntentType,'FOR PROCESSING')
                RawPro='RAW';
            elseif strcmp(metadata.PresentationIntentType,'FOR PRESENTATION')
                RawPro='PRO';
            else
                RawPro='--';
            end
        else
            RawPro='--';
        end
        
        if isfield(metadata,'PatientID')
            NHSNumber=metadata.PatientID;
        else
            NHSNumber='--';
        end
        if isfield(metadata,'PatientBirthDate')
            DoB=metadata.PatientBirthDate;
        else
            DoB='--';
        end
        if isfield(metadata,'PatientAge')
            Age=metadata.PatientAge;
        else
            Age='--';
        end
        if isfield(metadata,'AccessionNumber')
            AccNo=metadata.AccessionNumber;
        else
            AccNo='--';
        end
        
        if isfield(metadata,'PatientName')
           if isfield(metadata.PatientName,'GivenName')
               PgName=metadata.PatientName.GivenName;
           else
               PgName='--';
           end
           if isfield(metadata.PatientName,'FamilyName')
               PfName=metadata.PatientName.FamilyName;
           else
               PfName='--';
           end
        else
            PgName='--';
            PfName='--';
        end
        
        if isfield(metadata,'ReferringPhysicianName')
           if isfield(metadata.ReferringPhysicianName,'GivenName')
               RgName=metadata.ReferringPhysicianName.GivenName;
           else
               RgName='--';
           end
           if isfield(metadata.ReferringPhysicianName,'FamilyName')
               RfName=metadata.ReferringPhysicianName.FamilyName;
           else
               RfName='--';
           end
        else
            RgName='--';
            RfName='--';
        end
        
        if isfield(metadata,'PerformingPhysicianName')
           if isfield(metadata.PerformingPhysicianName,'FamilyName')
               PPfName=metadata.PerformingPhysicianName.FamilyName;
           else
               PPfName='--';
           end
        else
            PPfName='--';
        end
        
        if isfield(metadata,'OperatorName')
           if isfield(metadata.OperatorName,'FamilyName')
               OfName=metadata.OperatorName.FamilyName;
           else
               OfName='--';
           end
        else
            OfName='--';
        end
        % name the new DICOM file
        infoCell{numImage,1}=' ';
        infoCell{numImage,2}=' ';
        infoCell{numImage,3}=' ';
        infoCell{numImage,4}=' ';
        infoCell{numImage,5}=' ';
        infoCell{numImage,6}=' ';
        infoCell{numImage,7}=' ';
        infoCell{numImage,8}=' ';
        infoCell{numImage,9}=' ';
              
        infoCell{numImage,1}=NHSNumber;
        infoCell{numImage,2}=scanDate;
        infoCell{numImage,3}=DoB;
        infoCell{numImage,4}=Age;
        infoCell{numImage,5}=AccNo;
        infoCell{numImage,6}=[PgName,'_', PfName];
        infoCell{numImage,7}=[RgName,'_', RfName];
        infoCell{numImage,8}=[PPfName,'/',OfName];
        infoCell{numImage,9}=[databaseName,'-',subStr,'-',LR,View,sTime,'-', scanDate,'-',RawPro,'.dcm'];
                
        % Anonoymous confidential information
        if isfield(metadata,'PatientID')
            metadata.PatientID=[databaseName,'-',subStr];
        end
        if isfield(metadata,'AccessionNumber')
            metadata.AccessionNumber='';
        end
        if isfield(metadata,'PatientName')
            if isfield(metadata.PatientName,'FamilyName')
                metadata.PatientName.FamilyName='';
            end
        end
        if isfield(metadata,'PatientName')
            if isfield(metadata.PatientName,'GivenName')
                metadata.PatientName.GivenName='';
            end
        end
        if isfield(metadata,'PatientName')
            if isfield(metadata.PatientName,'MiddleName')
                metadata.PatientName.MiddleName='';
            end
        end
        if isfield(metadata,'PatientName')
            if isfield(metadata.PatientName,'NamePrefix')
                metadata.PatientName.NamePrefix='';
            end
        end
        if isfield(metadata,'PatientName')
            if isfield(metadata.PatientName,'NameSuffix')
                metadata.PatientName.NameSuffix='';
            end
        end
        if isfield(metadata,'PatientBirthDate')
            metadata.PatientBirthDate='';
        end
        if isfield(metadata,'ReferringPhysicianName')
            if isfield(metadata.ReferringPhysicianName,'FamilyName')
                metadata.ReferringPhysicianName.FamilyName='';
            end
        end
        if isfield(metadata,'ReferringPhysicianName')
            if isfield(metadata.ReferringPhysicianName,'GivenName')
                metadata.ReferringPhysicianName.GivenName='';
            end
        end
        if isfield(metadata,'ReferringPhysicianName')
            if isfield(metadata.ReferringPhysicianName,'MiddleName')
                metadata.ReferringPhysicianName.MiddleName='';
            end
        end
        if isfield(metadata,'ReferringPhysicianName')
            if isfield(metadata.ReferringPhysicianName,'NamePrefix')
                metadata.ReferringPhysicianName.NamePrefix='';
            end
        end
        if isfield(metadata,'ReferringPhysicianName')
            if isfield(metadata.ReferringPhysicianName,'NameSuffix')
                metadata.ReferringPhysicianName.NameSuffix='';
            end
        end
        if isfield(metadata,'PerformingPhysicianName')
            metadata.PerformingPhysicianName='';
        end
        if isfield(metadata,'OperatorName')
            metadata.OperatorName='';
        end
        if isfield(metadata,'BitsAllocated')
            metadata.BitsAllocated='14';
        end
        if isfield(metadata,'BitsStored')
            metadata.BitsStored='14';
        end
        if isfield(metadata,'HighBit')
            metadata.HighBit='13';
        end
        if isfield(metadata,'BodyPartExamined')
            metadata.BodyPartExamined='BREAST';
        end

        indxP=find(savePathname=='\');
        if ~exist([savePathname(1:indxP(2)),databaseName,'_',subStr])
            mkdir([savePathname(1:indxP(2)),databaseName,'_',subStr]);
        end
       % dicomanon(imageName, [savePathname(1:indxP(2)),databaseName,'_',subStr,'\', databaseName,'-',subStr,'-',LR,View,sTime,'-', scanDate,'-',RawPro,'.dcm'], 'keep','PatientAge');
        dicomwrite(X, [savePathname(1:indxP(2)),databaseName,'_',subStr,'\', databaseName,'-',subStr,'-',LR,View,sTime,'-', scanDate,'-LND',NHSNumber(end-2:end),'-',RawPro,'.dcm'], metadata,'CreateMode', 'copy'); %
      % if exist(savePathname)~=0 %&& length(find(savePathname=='\'))>2,
       %     rmdir(savePathname);
       % end
        
        numImage=numImage+1;
        
        pastImnum=pastImnum+1;
        waitbar( pastImnum/totImage , w);
    end
end
e = etime(clock, t);
%subjectID=subjectID+1;

end


