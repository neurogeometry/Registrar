function varargout = Registrar(varargin)
% ============================== About ====================================
% ----------------Copyright 2019 Northeastern University-------------------
%
% Purpose: Registrar GUI handler
% -------------------------------------------------------------------------
% Author: Seyed Mostafa Mousavi Kahaki and Armen Stepanyants
% Northeastern University, Boston, MA, USA
% kahaki@neu.edu, a.stepanyants.neu.edu
% =========================================================================
% -------------------------------------------------------------------------
% Registrar MATLAB code for Registrar.fig
%      Registrar, by itself, creates a new Registrar or
%      raises the existing singleton*.
%
%      H = Registrar returns the handle to a new Registrar or
%      the handle to the existing singleton*.
%
%      Registrar('CALLBACK',hObject,eventData,handles,...) calls the
%      local function named CALLBACK in Registrar.M with the given
%      input arguments.
%
%      Registrar('Property','Value',...) creates a new
%      Registrar or raises the existing singleton*.  Starting from
%      the left, propertyg value pairs are applied to the GUI before
%      Registrar_OpeningFcn gets called.  An unrecognized property
%      name or invalid value makes property application stop.  All inputs
%      are passed to Registrar_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Registrar

% Last Modified by GUIDE v2.5 01-Feb-2019 15:17:37

% Begin initialization code - DO NOT EDIT
clc;
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Registrar_OpeningFcn, ...
    'gui_OutputFcn',  @Registrar_OutputFcn, ...
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


% --- Executes just before Registrar is made visible.
function Registrar_OpeningFcn(hObject, eventdata, handles, varargin)
% Update Icon
try
    jFrame=get(handles.figure1,'javaframe');
    jicon=javax.swing.ImageIcon('icon.png');
    jFrame.setFigureIcon(jicon);
catch
end

handles.output = hObject;
guidata(hObject, handles);
handles.output = hObject;
myCluster = parcluster('local');
MaxNumWorkers = myCluster.NumWorkers;
set(handles.text28,'Enable','off');
set(handles.text28, 'String', ['(max: ',num2str(MaxNumWorkers),')']);
set(handles.edit14, 'String', num2str(MaxNumWorkers));
% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = Registrar_OutputFcn(hObject, eventdata, handles)
if sum(size(handles.axes3.Children)) == 0
    patch('XData',[0,0,100,100],'YData',[0,20,20,0],'FaceColor','white','Parent',handles.axes3);
end
% Get default command line output from handles structure
varargout{1} = handles.output;

function edit1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pushbutton7_Callback(hObject, eventdata, handles)
addpath('Functions');
parameters;
TransformationValue = get(handles.popupmenu3,'Value'); % 1=Translation | 2=Rigis | 3=Affine | 4=NonRigid
Seq_Par = get(handles.popupmenu4,'Value'); % 1=Sequential | 2=Parallel
Par_workers = str2double(get(handles.edit14,'String')); % Number of Workers
blendingSID = str2double(get(handles.edit15,'String'));
StackList_csv_pth = get(handles.edt_stacklist,'String');
if handles.checkbox15.Value
    Log();
    LogHandle=findobj(0,'Name','Log');
    LogHandle.Children(2).String = {};
else
    LogHandle = [];
end

tic
registeration(StackList_csv_pth,TransformationValue,Seq_Par,Par_workers,blendingSID,handles,LogHandle)
% try
%     registeration(StackList_csv_pth,TransformationValue,Seq_Par,Par_workers,blendingSID,handles,LogHandle)
% catch ME
%     LogHandle.Children(2).String = ME.getReport;
% end
TotalTime = toc

function pushbutton1_Callback(hObject, eventdata, handles)
[filename, pathname] = uigetfile('*.csv', 'Pick a Stack List CSV file');
if isequal(filename,0) || isequal(pathname,0)
    disp('User pressed cancel')
else
    disp(['User selected ', fullfile(pathname, filename)])
    set(handles.edt_stacklist, 'String', fullfile(pathname, filename));
end

function pushbutton4_Callback(hObject, eventdata, handles)
[filename, pathname] = uigetfile('*.csv', 'Pick a Stack Positions CSV file');
if isequal(filename,0) || isequal(pathname,0)
    disp('User pressed cancel')
else
    disp(['User selected ', fullfile(pathname, filename)])
    set(handles.edit9, 'String', fullfile(pathname, filename));
    StackPositions_pixels = xlsread(fullfile(pathname, filename));
    set(handles.listbox3, 'String', num2str(StackPositions_pixels));
end

function pushbutton5_Callback(hObject, eventdata, handles)
[filename, pathname] = uigetfile('*.csv', 'Pick a Stack Overlap CSV file');
if isequal(filename,0) || isequal(pathname,0)
    disp('User pressed cancel')
else
    disp(['User selected ', fullfile(pathname, filename)])
    set(handles.edit10, 'String', fullfile(pathname, filename));
end

function pushbutton6_Callback(hObject, eventdata, handles)
[filename, pathname] = uigetfile('*.csv', 'Pick a Stack Sizes CSV file');
if isequal(filename,0) || isequal(pathname,0)
    disp('User pressed cancel')
else
    disp(['User selected ', fullfile(pathname, filename)])
    set(handles.edit11, 'String', fullfile(pathname, filename));
end

function pushbutton8_Callback(hObject, eventdata, handles)
winopen('data/StackPositions_Registered.csv');

function pushbutton9_Callback(hObject, eventdata, handles)
winopen('data/StackPositions_Registered.csv');

% --------------------------------------------------------------------
function docmnu_Callback(hObject, eventdata, handles)
docurl = 'http://neurogeometry.com/';
web(docurl);

% --------------------------------------------------------------------
function websitemnu_Callback(hObject, eventdata, handles)
docurl = 'http://neurogeometry.com/';
web(docurl);

function v_Callback(hObject, eventdata, handles)
val = get(handles.v,'Value');

if val == 1 %MultiStack Registration
    set(handles.edt_stacklist,'String','../Data/Neocortical.csv');
    handles.chkretilling.Enable = 'on';
elseif val == 2 % TimeLapse Registration
    set(handles.edt_stacklist,'String','../Data/DL083.csv');
    handles.chkretilling.Value = 0;
    handles.chkretilling.Enable = 'off';
elseif val == 3 % Stack Registration
    set(handles.edt_stacklist,'String','../Data/DL083B001G.csv');
    handles.chkretilling.Value = 0;
    handles.chkretilling.Enable = 'off';
end

function v_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function checkbox6_Callback(hObject, eventdata, handles)
val = get(handles.checkbox6,'Value');
if val == 0
    set(handles.checkbox13,'Enable','off');
    set(handles.checkbox14,'Enable','off');
    set(handles.checkbox10,'Enable','off');
else
    set(handles.checkbox13,'Enable','on');
    set(handles.checkbox14,'Enable','on');
    set(handles.checkbox10,'Enable','on');
end

function popupmenu4_Callback(hObject, eventdata, handles)
val = get(handles.popupmenu4,'Value');
if val == 1
    set(handles.edit14,'Enable','off');
    set(handles.text27,'Visible','off');
    set(handles.text26,'Enable','off');
    set(handles.text28,'Enable','off');
else
    set(handles.edit14,'Enable','on');
    set(handles.text27,'Visible','on');
    set(handles.text26,'Enable','on');
    set(handles.text28,'Enable','on');
end

function edit14_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit15_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function checkbox8_Callback(hObject, eventdata, handles)
val = get(handles.checkbox8,'Value');
if val
    Tile3D1=evalin('base','Tile3D');
    volumeViewer(Tile3D1);
else
    volumeViewer close
end

function checkbox10_Callback(hObject, eventdata, handles)
val = get(handles.checkbox10,'Value');
if val
    try
        Tile3D1=evalin('base','Tile3D');
        volumeViewer(Tile3D1);
    catch
        warndlg('Can not load tiles, Please run Global Optimization and check Preview','!! Warning !!');
        volumeViewer CLOSE
    end
else
    volumeViewer CLOSE
end

function pushbutton10_Callback(hObject, eventdata, handles)
set(handles.pushbutton10,'userdata',1)

disp('Stop');
return

function pushbutton10_KeyPressFcn(hObject, eventdata, handles)
set(handles.pushbutton10,'userdata',1)
disp('Stop');
return

function chkdebug_Callback(hObject, eventdata, handles)

if handles.chkdebug.Value
    Debug();
    DebugHandle=findobj(0,'Name','Debug');
    DebugHandle.Visible = 'on';
else
    try
        Debug();
        DebugHandle=findobj(0,'Name','Debug');
        DebugHandle.Visible = 'off';
    catch
    end
end

function chkretilling_Callback(hObject, eventdata, handles)
val = get(handles.chkretilling,'Value');
if val == 1
    set(handles.popupretilling,'Enable','on');
    set(handles.txtoutput,'Enable','on');
else
    set(handles.popupretilling,'Enable','off');
    set(handles.txtoutput,'Enable','off');
end

% --------------------------------------------------------------------
function parameters_mnu_Callback(hObject, eventdata, handles)
!notepad Functions\parameters.m
% open('Functions\parameters.m');

% --------------------------------------------------------------------
function website_mnu_Callback(hObject, eventdata, handles)
web('http://www.northeastern.edu/neurogeometry/resources/registrar');

% --------------------------------------------------------------------
function manual_mnu_Callback(hObject, eventdata, handles)
web('http://www.northeastern.edu/neurogeometry/resources/registrar');
%open('UserManual.docx');

function checkbox13_Callback(hObject, eventdata, handles)
if handles.checkbox13.Value
    Visualization();
    VisualizationHandle=findobj(0,'Name','Visualization');
    Axes1V = VisualizationHandle.Children(3);
    ButtonGroup1V = VisualizationHandle.Children(2);
    AfterButton = ButtonGroup1V.Children(1);
    set(ButtonGroup1V,'SelectedObject',AfterButton);
    try
        Tile3D=evalin('base','Tile3D');
        h_im=imshow(max(Tile3D,[],3),[0 max(Tile3D(:))],'Parent',Axes1V);
        VisualizationHandle.Children(3).Children.CDataMapping = 'direct';
        VisualizationHandle.Visible = 'on';
    catch
        warndlg('Can not load tiles, Please run Global Optimization and check Preview','!! Warning !!');
        VisualizationHandle.Visible = 'off';
    end
else
    try
        VisualizationHandle=findobj(0,'Name','Visualization');
        VisualizationHandle.Visible = 'off';
    catch
    end
end

function checkbox14_Callback(hObject, eventdata, handles)

if handles.checkbox14.Value
    VisualiztionStack();
    VisualizationStackHandle=findobj(0,'Name','VisualiztionStack');
    Axes1V = VisualizationStackHandle.Children(3);
    ButtonGroup1V = VisualizationStackHandle.Children(2);
    AfterButton = ButtonGroup1V.Children(1);
    set(ButtonGroup1V,'SelectedObject',AfterButton);
    try
        Tile3D=evalin('base','Tile3D');
        plainsViewer(VisualizationStackHandle,Tile3D,1);
        VisualizationStackHandle.Visible = 'on';
        VisualizationStackHandle.Children(3).Children.CDataMapping = 'direct';
    catch
        warndlg('Can not load tiles, Please run Global Optimization and check Preview','!! Warning !!');
        VisualizationStackHandle.Visible = 'off';
    end
else
    try
        VisualizationStackHandle=findobj(0,'Name','VisualiztionStack');
        VisualizationStackHandle.Visible = 'off';
    catch
    end
end

function figure1_DeleteFcn(hObject, eventdata, handles)
delete(findall(0));

function figure1_CloseRequestFcn(hObject, eventdata, handles)
delete(hObject);

function checkbox15_Callback(hObject, eventdata, handles)
if handles.checkbox15.Value
    Log();
    LogHandle=findobj(0,'Name','Log');
    LogHandle.Visible = 'on';
else
    try
        Log();
        LogHandle=findobj(0,'Name','Log');
        LogHandle.Visible = 'off';
    catch
    end
end

function checkbox16_Callback(hObject, eventdata, handles)
if handles.checkbox16.Value
    DatasetMap();
    MapHandle=findobj(0,'Name','Stack Map');
    MapHandle.Visible = 'on';
else
    try
        DatasetMap();
        MapHandle=findobj(0,'Name','Stack Map');
        MapHandle.Visible = 'off';
    catch
    end
end
