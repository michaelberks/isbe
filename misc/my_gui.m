function varargout = my_gui(varargin)
% MY_GUI M-file for my_gui.fig
%      MY_GUI, by itself, creates a new MY_GUI or raises the existing
%      singleton*.
%
%      H = MY_GUI returns the handle to a new MY_GUI or the handle to
%      the existing singleton*.
%
%      MY_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MY_GUI.M with the given input arguments.
%
%      MY_GUI('Property','Value',...) creates a new MY_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before my_gui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to my_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help my_gui

% Last Modified by GUIDE v2.5 16-Aug-2006 13:00:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @my_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @my_gui_OutputFcn, ...
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


% --- Executes just before my_gui is made visible.
function my_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to my_gui (see VARARGIN)

% Choose default command line output for my_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
handles.mass_border = [];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes my_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = my_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --------------------------------------------------------------------
function file_menu_Callback(hObject, eventdata, handles)
% hObject    handle to file_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function menu_file_open_orig_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_open_orig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename pathname] = uigetfile('*.bmp','File Selector');
handles.im_orig = imread(strcat(pathname, filename));
set(handles.file_orig_text, 'String', strcat(pathname, filename));

handles.fig1 = openfig('annotate1'); colormap('gray');
handles.axes1 = axes;
set(handles.fig1,'CurrentAxes', handles.axes1)
axes(handles.axes1)
handles.im1 = imagesc(handles.im_orig); axis image; hold on;
set(handles.im1, 'ButtonDownFcn', @plot_point_Callback1);

handles.fig2 = openfig('annotate2'); colormap('jet');
handles.axes2 = axes;
set(handles.fig2,'CurrentAxes', handles.axes2)
axes(handles.axes2)
handles.im2 = imagesc(handles.im_orig); axis image; hold on;
set(handles.im2, 'ButtonDownFcn', @plot_point_Callback2);

handles.fig3 = openfig('annotate3'); colormap('spring');
handles.axes3 = axes;
set(handles.fig3,'CurrentAxes', handles.axes3)
axes(handles.axes3)
handles.im3 = imagesc(handles.im_orig); axis image; hold on;
set(handles.im3, 'ButtonDownFcn', @plot_point_Callback3);

clear handles.im_orig;

% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function menu_file_save_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_edit_copy_Callback(hObject, eventdata, handles)
% hObject    handle to menu_edit_copy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_edit_cut_Callback(hObject, eventdata, handles)
% hObject    handle to menu_edit_cut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_edit_paste_Callback(hObject, eventdata, handles)
% hObject    handle to menu_edit_paste (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_edit_Callback(hObject, eventdata, handles)
% hObject    handle to menu_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function file_orig_text_Callback(hObject, eventdata, handles)
% hObject    handle to file_orig_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of file_orig_text as text
%        str2double(get(hObject,'String')) returns contents of file_orig_text as a double


% --- Executes during object creation, after setting all properties.
function file_orig_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to file_orig_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --------------------------------------------------------------------
function plot_point_Callback1(hObject, eventdata)
% hObject    handle to plot_point_Callback (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles
XYZ = get(handles.axes1, 'CurrentPoint');
handles.mass_border = [handles.mass_border; [XYZ(1,1), XYZ(1,2)]];

plot_point(handles, XYZ(1,1), XYZ(1,2));
display('This worked!');

% --------------------------------------------------------------------
function plot_point_Callback2(hObject, eventdata)
% hObject    handle to plot_point_Callback (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles
XYZ = get(handles.axes2, 'CurrentPoint');
handles.mass_border = [handles.mass_border; [XYZ(1,1), XYZ(1,2)]];

plot_point(handles, XYZ(1,1), XYZ(1,2));

% Update handles structure
guidata(hObject, handles);
display('This worked!');

% --------------------------------------------------------------------
function plot_point_Callback3(hObject, eventdata)
% hObject    handle to plot_point_Callback (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles
XYZ = get(handles.axes3, 'CurrentPoint');
handles.mass_border = [handles.mass_border; [XYZ(1,1), XYZ(1,2)]];

plot_point(handles, XYZ(1,1), XYZ(1,2));

% Update handles structure
guidata(hObject, handles);
display('This worked!');

% --------------------------------------------------------------------
function plot_point(handles, x, y)

figure(handles.fig1);
axes(handles.axes1);
plot(x, y, 'rx');

figure(handles.fig2);
axes(handles.axes2);
plot(x, y, 'rx');

figure(handles.fig3);
axes(handles.axes3);
plot(x, y, 'rx');


