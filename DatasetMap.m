function varargout = DatasetMap(varargin)
% DATASETMAP MATLAB code for DatasetMap.fig
%      DATASETMAP, by itself, creates a new DATASETMAP or raises the existing
%      singleton*.
%
%      H = DATASETMAP returns the handle to a new DATASETMAP or the handle to
%      the existing singleton*.
%
%      DATASETMAP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DATASETMAP.M with the given input arguments.
%
%      DATASETMAP('Property','Value',...) creates a new DATASETMAP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DatasetMap_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DatasetMap_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DatasetMap

% Last Modified by GUIDE v2.5 25-Jan-2019 15:06:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DatasetMap_OpeningFcn, ...
                   'gui_OutputFcn',  @DatasetMap_OutputFcn, ...
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


% --- Executes just before DatasetMap is made visible.
function DatasetMap_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DatasetMap (see VARARGIN)
try
    jFrame=get(handles.figure1,'javaframe');
    jicon=javax.swing.ImageIcon('icon.png');
    jFrame.setFigureIcon(jicon);
catch
end
% Choose default command line output for DatasetMap
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DatasetMap wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DatasetMap_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
NCT_RegistrationHandle=findobj(0,'Name','Registrar');
NCT_RegistrationHandle.Children.findobj('Tag','checkbox16').Value = 0;
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
