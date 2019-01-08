function varargout = VisualiztionStack(varargin)
% VISUALIZTIONSTACK MATLAB code for VisualiztionStack.fig
%      VISUALIZTIONSTACK, by itself, creates a new VISUALIZTIONSTACK or raises the existing
%      singleton*.
%
%      H = VISUALIZTIONSTACK returns the handle to a new VISUALIZTIONSTACK or the handle to
%      the existing singleton*.
%
%      VISUALIZTIONSTACK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VISUALIZTIONSTACK.M with the given input arguments.
%
%      VISUALIZTIONSTACK('Property','Value',...) creates a new VISUALIZTIONSTACK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before VisualiztionStack_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to VisualiztionStack_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help VisualiztionStack

% Last Modified by GUIDE v2.5 08-Jan-2019 11:45:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @VisualiztionStack_OpeningFcn, ...
                   'gui_OutputFcn',  @VisualiztionStack_OutputFcn, ...
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


% --- Executes just before VisualiztionStack is made visible.
function VisualiztionStack_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to VisualiztionStack (see VARARGIN)

% Choose default command line output for VisualiztionStack
handles.output = hObject;
jFrame=get(handles.figure1,'javaframe');
jicon=javax.swing.ImageIcon('icon.png');
jFrame.setFigureIcon(jicon);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes VisualiztionStack wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = VisualiztionStack_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
NCT_RegistrationHandle=findobj(0,'Name','Registrar');
NCT_RegistrationHandle.Children(2).Children(9).Children(1).Value = 0;
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when selected object is changed in uibuttongroup1.
function uibuttongroup1_SelectionChangedFcn(hObject, eventdata, handles)
valBeforeAfter = get(handles.uibuttongroup1.SelectedObject,'Tag');

if strcmp(valBeforeAfter,'radiobutton1')
    Tile3D=evalin('base','Tile3D_org');

else
    Tile3D=evalin('base','Tile3D');

end
    VisualizationStackHandle=findobj(0,'Name','VisualiztionStack');
    set(VisualizationStackHandle.Children(1),'Visible','on');
    set(VisualizationStackHandle.Children(1),'Min',1,'Max',size(Tile3D,3),'Value', 1, 'SliderStep', [1/size(Tile3D,3) 1/size(Tile3D,3)]);
    %     set(handles.slider2,'Max',size(Tile3D1,3));
    %     set(handles.slider2,'SliderStep',1);
    plainsViewer(VisualizationStackHandle,Tile3D);
%     h_im=imshow(max(Tile3D,[],3),[0 max(Tile3D(:))],'Parent',Axes1V);
    VisualizationStackHandle.Visible = 'on';
% hObject    handle to the selected object in uibuttongroup1 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
