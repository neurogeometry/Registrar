function varargout = Debug(varargin)
% DEBUG MATLAB code for Debug.fig
%      DEBUG, by itself, creates a new DEBUG or raises the existing
%      singleton*.
%
%      H = DEBUG returns the handle to a new DEBUG or the handle to
%      the existing singleton*.
%
%      DEBUG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DEBUG.M with the given input arguments.
%
%      DEBUG('Property','Value',...) creates a new DEBUG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Debug_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Debug_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Debug

% Last Modified by GUIDE v2.5 07-Jan-2019 14:37:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Debug_OpeningFcn, ...
                   'gui_OutputFcn',  @Debug_OutputFcn, ...
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


% --- Executes just before Debug is made visible.
function Debug_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Debug (see VARARGIN)

% Choose default command line output for Debug
handles.output = hObject;
try
    jFrame=get(handles.figure1,'javaframe');
    jicon=javax.swing.ImageIcon('icon.png');
    jFrame.setFigureIcon(jicon);
catch
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Debug wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Debug_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
NCT_RegistrationHandle=findobj(0,'Name','Registrar');
NCT_RegistrationHandle.findobj('Tag','chkdebug').Value = 0;
% NCT_RegistrationHandle
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
