function varargout = GUIAnonymise(varargin)
% GUIANONYMISE M-file for GUIAnonymise.fig
%      GUIANONYMISE, by itself, creates a new GUIANONYMISE or raises the existing
%      singleton*.
%
%      H = GUIANONYMISE returns the handle to a new GUIANONYMISE or the handle to
%      the existing singleton*.
%
%      GUIANONYMISE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUIANONYMISE.M with the given input arguments.
%
%      GUIANONYMISE('Property','Value',...) creates a new GUIANONYMISE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUIAnonymise_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUIAnonymise_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUIAnonymise

% Last Modified by GUIDE v2.5 02-May-2013 14:27:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUIAnonymise_OpeningFcn, ...
                   'gui_OutputFcn',  @GUIAnonymise_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before GUIAnonymise is made visible.
function GUIAnonymise_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUIAnonymise (see VARARGIN)

% Choose default command line output for GUIAnonymise
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUIAnonymise wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUIAnonymise_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in AnonyDICOM.
function AnonyDICOM_Callback(hObject, eventdata, handles)
% hObject    handle to AnonyDICOM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

startSubjectID=str2double(get(handles.subjectID,'string'));
databaseName=get(handles.databaseName,'string');

dirName = uigetdir('C:\', 'Choose the image folder');
if (~ischar(dirName))
    return;
end

outPathname=cell(1,99999);
inPathname=cell(1,99999);

outUpdatePath = uigetdir('C:\', 'Choose an output folder');
inUpdatePath=dirName;
outPathname{1}=outUpdatePath;
inPathname{1}=inUpdatePath;
nFolder=2;

w = waitbar(0, 'In progress...');

AllFileName = ls ( sprintf('%s\\',inUpdatePath) );
for i=1:size(AllFileName,1)
    indF=(AllFileName(i,:)~=' ');
    tempInd=find(indF==1);
    indF(1:max(tempInd))=1;
    B=AllFileName(i,indF);
    if ~strcmp(B,'.') && ~strcmp(B,'..') && isdir([inUpdatePath,'\',B]),
       inPathname{nFolder}=[inUpdatePath,'\',B];
       outPathname{nFolder}=[outUpdatePath,'\',B]; 
    %   mkdir(outUpdatePath,B);
       nFolder=nFolder+1;
       tmpIn=[inUpdatePath,'\',B];
       tmpOut=[outUpdatePath,'\',B]; 
       [inPathname, outPathname, nFolder]=checkAllFolder(inPathname, outPathname, nFolder,tmpIn,tmpOut);
    end
    
end

totImage=0;
for i=1:(nFolder-1)
    inFolder=inPathname{i};
    imageFilename = ls ( sprintf('%s\\*.*',inFolder) );
    totImage = totImage+ size(imageFilename,1);
end

infoCell=cell(totImage,9);
infoCell{1,1}='PatientID';
infoCell{1,2}='StudyDate';
infoCell{1,3}='PatientBirthDate';
infoCell{1,4}='PatientAge';
infoCell{1,5}='AccessionNumber';
infoCell{1,6}='PatientName';
infoCell{1,7}='PhysicianName';
infoCell{1,8}='PerformingPhysicianName/Operator';
infoCell{1,9}='ImageName';
    
totTime=0;
pastImnum=0;
subjectID=startSubjectID;
numImage=2;

inUpdatePath=dirName;

for i=1:size(AllFileName,1)
    indF=(AllFileName(i,:)~=' ');
    tempInd=find(indF==1);
    indF(1:max(tempInd))=1;
    B=AllFileName(i,indF);
    outPathname=cell(1,99999);
    inPathname=cell(1,99999);
    outPathname{1}=outUpdatePath;
    inPathname{1}=inUpdatePath;
    
    if ~strcmp(B,'.') && ~strcmp(B,'..') && isdir([inUpdatePath,'\',B]),
       nFolder=2;
       inPathname{nFolder}=[inUpdatePath,'\',B];
       outPathname{nFolder}=[outUpdatePath,'\',B]; 
    %   mkdir(outUpdatePath,B);
       nFolder=nFolder+1;
       tmpIn=[inUpdatePath,'\',B];
       tmpOut=[outUpdatePath,'\',B];
       [inPathname, outPathname, nFolder]=checkAllFolder(inPathname, outPathname, nFolder,tmpIn,tmpOut);
       for j=1:(nFolder-1)
           inFolder=inPathname{j};
           [proc_time, pastImnum, infoCell,numImage] = AnonmousDICOM (inFolder, [outPathname{j},'\'], totImage, pastImnum, w, subjectID, databaseName, infoCell,numImage);
           totTime=totTime+proc_time;
           waitbar( pastImnum/totImage , w);
       end
       subjectID=subjectID+1;
    end
    
end

close(w);

xlswrite([outUpdatePath,'\',databaseName,'_confidential.xls'], infoCell);

resMessage = sprintf ('The final subject ID is %8d. Elapsed time to analyse all images is %8.2f seconds.'...
    ,subjectID-1, totTime);
helpdlg(resMessage,'Analysis completed');


guidata(hObject, handles);



function databaseName_Callback(hObject, eventdata, handles)
% hObject    handle to databaseName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of databaseName as text
%        str2double(get(hObject,'String')) returns contents of databaseName as a double


% --- Executes during object creation, after setting all properties.
function databaseName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to databaseName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function subjectID_Callback(hObject, eventdata, handles)
% hObject    handle to subjectID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of subjectID as text
%        str2double(get(hObject,'String')) returns contents of subjectID as a double


% --- Executes during object creation, after setting all properties.
function subjectID_CreateFcn(hObject, eventdata, handles)
% hObject    handle to subjectID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
