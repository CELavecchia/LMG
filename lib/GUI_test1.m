function [varargout,model] = GUI_test1(varargin)
% GUI_TEST1 MATLAB code for GUI_test1.fig
%      GUI_TEST1, by itself, creates a new GUI_TEST1 or raises the existing
%      singleton*.
%
%      H = GUI_TEST1 returns the handle to a new GUI_TEST1 or the handle to
%      the existing singleton*.
%
%      GUI_TEST1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_TEST1.M with the given input arguments.
%
%      GUI_TEST1('Property','Value',...) creates a new GUI_TEST1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_test1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_test1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_test1

% Last Modified by GUIDE v2.5 15-Jun-2017 19:54:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_test1_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_test1_OutputFcn, ...
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


% --- Executes just before GUI_test1 is made visible.
function GUI_test1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_test1 (see VARARGIN)

% Choose default command line output for GUI_test1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_test1 wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_test1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles;
handles.output = hObject;

%function GUI_test1_OutputFnc(hObject, eventdata, handles) 
%guidata(hObject, handles)

% Update handles structure
%guidata(hObject, handles);


% --- Executes on button press in radiobutton1.
function model = radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if hObject.Value ==1
    model.type = 'average';
end
%handles.output = hObject;

% Update handles structure
%guidata(hObject, handles);
    % Hint: get(hObject,'Value') returns toggle state of radiobutton1
%varargout{2} = get(hObject,'Value');

% --- Executes on button press in radiobutton2.
function model =radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if hObject.Value ==1
    model.type = 'subject_specific';
end
%handles.output = hObject;

% Update handles structure
%guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of radiobutton2
%varargout{3} = get(hObject,'Value');


function model =edit1_Callback(hObject, eventdata, handles)%age
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% model.average.age = hObject.Value;
% handles.output = hObject;

%guidata(hObject, handles);


% Update handles structure
%guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
%varargout{4} = str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles) %height
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
model.average.height = hObject.Value;
handles.output = hObject;

% Update handles structure
%guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
%varargout{5} = str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function model =edit4_Callback(hObject, eventdata, handles) %lumbar curvature
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if hObject.Value==0
    model.average.angle = 43;%degree
else 
    model.average.angle = hObject.Value;
end
    

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double
%varargout{6} = str2double(get(hObject,'String'))

% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton3.
function model =radiobutton3_Callback(hObject, eventdata, handles) %male
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if hObject.Value == 1;
    model.average.sex =1;
end
%handles.output = hObject;

% Update handles structure
%guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of radiobutton3
%varargout{7} = get(hObject,'Value');

% --- Executes on button press in radiobutton4.
function model =radiobutton4_Callback(hObject, eventdata, handles) %female
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if hObject.Value == 1;
    model.average.sex = 0;
end
handles.output = hObject;

% Update handles structure
%guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of radiobutton4
%varargout{8} = get(hObject,'Value');

%varargout{1} = hObject;
%varargout{2} = handles.string;
%save 'guioutput'
