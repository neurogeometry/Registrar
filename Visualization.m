function varargout = Visualization(varargin)
% VISUALIZATION MATLAB code for Visualization.fig
%      VISUALIZATION, by itself, creates a new VISUALIZATION or raises the existing
%      singleton*.
%
%      H = VISUALIZATION returns the handle to a new VISUALIZATION or the handle to
%      the existing singleton*.
%
%      VISUALIZATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VISUALIZATION.M with the given input arguments.
%
%      VISUALIZATION('Property','Value',...) creates a new VISUALIZATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Visualization_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Visualization_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Visualization

% Last Modified by GUIDE v2.5 07-Jan-2019 15:01:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Visualization_OpeningFcn, ...
                   'gui_OutputFcn',  @Visualization_OutputFcn, ...
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


% --- Executes just before Visualization is made visible.
function Visualization_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Visualization (see VARARGIN)

% Choose default command line output for Visualization
handles.output = hObject;
jFrame=get(handles.figure1,'javaframe');
jicon=javax.swing.ImageIcon('icon.png');
jFrame.setFigureIcon(jicon);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Visualization wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Visualization_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes when selected object is changed in uibuttongroup1.
function uibuttongroup1_SelectionChangedFcn(hObject, eventdata, handles)
valBeforeAfter = get(handles.uibuttongroup1.SelectedObject,'Tag');

if strcmp(valBeforeAfter,'Before')
    Tile=evalin('base','Tile3D_org');
    
        P_XLim = handles.axes1.XLim;
        P_YLim = handles.axes1.YLim;
        P_Position = handles.axes1.Position;
        h_im=imshow(max(Tile,[],3),[0 max(Tile(:))],'Parent',handles.axes1);
        handles.axes1.XLim = P_XLim;
        handles.axes1.YLim = P_YLim;
        handles.axes1.Position=  P_Position;
    
else
    Tile=evalin('base','Tile3D');
    
        P_XLim = handles.axes1.XLim;
        P_YLim = handles.axes1.YLim;
        P_Position = handles.axes1.Position;
        h_im=imshow(max(Tile,[],3),[0 max(Tile(:))],'Parent',handles.axes1);
        handles.axes1.XLim = P_XLim;
        handles.axes1.YLim = P_YLim;
        handles.axes1.Position=  P_Position;

%         plainsViewer(handles,Tile);
   
end
% hObject    handle to the selected object in uibuttongroup1 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% % --- Executes when selected object is changed in uibuttongroup2.
% function uibuttongroup2_SelectionChangedFcn(hObject, eventdata, handles)
% valBeforeAfter = get(handles.uibuttongroup1.SelectedObject,'Tag');
% 
% if strcmp(valBeforeAfter,'Before')
%     Tile=evalin('base','Tile3D_org');
%     
%         P_XLim = handles.axes1.XLim;
%         P_YLim = handles.axes1.YLim;
%         P_Position = handles.axes1.Position;
%         h_im=imshow(max(Tile,[],3),[0 max(Tile(:))],'Parent',handles.axes1);
%         handles.axes1.XLim = P_XLim;
%         handles.axes1.YLim = P_YLim;
%         handles.axes1.Position=  P_Position;
%     
% else
%     Tile=evalin('base','Tile3D');
%     
%         P_XLim = handles.axes1.XLim;
%         P_YLim = handles.axes1.YLim;
%         P_Position = handles.axes1.Position;
%         h_im=imshow(max(Tile,[],3),[0 max(Tile(:))],'Parent',handles.axes1);
%         handles.axes1.XLim = P_XLim;
%         handles.axes1.YLim = P_YLim;
%         handles.axes1.Position=  P_Position;
% 
% %         plainsViewer(handles,Tile);
%    
% end
% hObject    handle to the selected object in uibuttongroup2 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


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
NCT_RegistrationHandle.Children(2).Children(9).Children(2).Value = 0;
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
