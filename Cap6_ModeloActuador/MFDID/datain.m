function varargout = datain(varargin)
% DATAIN Application M-file for datain.fig
%   DATAIN, by itself, creates a new DATAIN or raises the existing
%   singleton*.
%
%   H = DATAIN returns the handle to a new DATAIN or the handle to
%   the existing singleton*.
%
%   DATAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%   function named CALLBACK in DATAIN.M with the given input arguments.
%
%   DATAIN('Property','Value',...) creates a new DATAIN or raises the
%   existing singleton*.  Starting from the left, property value pairs are
%   applied to the GUI before datain_OpeningFunction gets called.  An
%   unrecognized property name or invalid value makes property application
%   stop.  All inputs are passed to datain_OpeningFcn via varargin.
%
%   *See GUI Options - GUI allows only one instance to run (singleton).
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help datain

% Last Modified by GUIDE v2.5 13-Mar-2015 14:23:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',          mfilename, ...
                   'gui_Singleton',     gui_Singleton, ...
                   'gui_OpeningFcn',    @datain_OpeningFcn, ...
                   'gui_OutputFcn',     @datain_OutputFcn, ...
                   'gui_LayoutFcn',     [], ...
                   'gui_Callback',      []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    varargout{1:nargout} = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before datain is made visible.
function datain_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to datain (see VARARGIN)

% Choose default command line output for datain
handles.output = [];

if nargin > 3,
	handles.ifile = varargin{1};
	evalin('base',['load ',handles.ifile]);
end

% Determine the position of the dialog - centered on the callback figure
% if available, else, centered on the screen
FigPos=get(0,'DefaultFigurePosition');
OldUnits = get(hObject, 'Units');
set(hObject, 'Units', 'pixels');
OldPos = get(hObject,'Position');
FigWidth = OldPos(3);
FigHeight = OldPos(4);
if isempty(gcbf)
    ScreenUnits=get(0,'Units');
    set(0,'Units','pixels');
    ScreenSize=get(0,'ScreenSize');
    set(0,'Units',ScreenUnits);

    FigPos(1)=1/2*(ScreenSize(3)-FigWidth);
    FigPos(2)=2/3*(ScreenSize(4)-FigHeight);
else
    GCBFOldUnits = get(gcbf,'Units');
    set(gcbf,'Units','pixels');
    GCBFPos = get(gcbf,'Position');
    set(gcbf,'Units',GCBFOldUnits);
    FigPos(1:2) = [(GCBFPos(1) + GCBFPos(3) / 2) - FigWidth / 2, ...
                   (GCBFPos(2) + GCBFPos(4) / 2) - FigHeight / 2];
end
FigPos(3:4)=[FigWidth FigHeight];
set(hObject, 'Position', FigPos);
set(hObject, 'Units', OldUnits);

% Updates the listbox to match the current workspace
vars = evalin('base','who');
set(handles.listbox_idata,'String',vars)

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes datain wait for user response (see UIRESUME)
uiwait(handles.datain);


% --- Outputs from this function are returned to the command line.
function varargout = datain_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

delete(handles.datain);


function varargout = f_OK_button_Callback(hObject, eventdata, handles)
% hObject    handle to f_OK_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.output.ind_SR_o = ...
    str2double(get(handles.edit_ind_SR_o,'String'));
handles.output.ind_SR_i = ...
    str2double(get(handles.edit_ind_SR_i,'String'));

try,
	handles.output.freq = evalin('base',get(handles.text_freq,'String'));
	handles.output.H_exp = evalin('base',get(handles.text_H_exp,'String'));
	
	try
		Sig_H = evalin('base',get(handles.text_Sig_H,'String'));
		Sig_H = Sig_H;
	catch
		try,
			coh = evalin('base',get(handles.text_coh,'String'));
			coh = coh;
	    catch
			coh = 0.9*ones(size(handles.output.H_exp));
	    end
	end
	try
		Sig_H_inv = evalin('base',get(handles.text_Sig_H_inv,'String'));
		Sig_H_inv = Sig_H_inv;
	end
	
	try,
	    handles.output.W_H = (1./Sig_H).^(1/2);		clear Sig_H
	catch,
	    handles.output.W_H = 1 ./ ...
	    	(sqrt(1-coh.^2) ./ abs(coh) .* abs(handles.output.H_exp));
	end
	if 0,
	W_max = max(max(handles.output.W_H(isfinite(handles.output.W_H))));
	handles.output.W_H(~isfinite(handles.output.W_H)) = W_max;
	end
	handles.output.W_H(~isfinite(handles.output.W_H)) = 0;
	
	try,
	    handles.output.W_H_inv = (1./Sig_H_inv).^(1/2);	clear Sig_H_inv
	catch,
	    handles.output.W_H_inv = 1 ./ ...
	    	(sqrt(1-coh.^2) ./ abs(coh) ./ abs(handles.output.H_exp));
	end
	clear coh
	if 0,
	W_inv_max = ...
		max(max(handles.output.W_H_inv(isfinite(handles.output.W_H_inv))));
	handles.output.W_H_inv(~isfinite(handles.output.W_H_inv)) = W_inv_max;
	end
	handles.output.W_H_inv(~isfinite(handles.output.W_H_inv)) = 0;
end

handles.output.n_o = ...
    str2double(get(handles.edit_n_o,'String'));
handles.output.n_i = ...
    str2double(get(handles.edit_n_i,'String'));

handles.output.n_pole = ...
    str2double(get(handles.edit_n_pole,'String'));
handles.output.n_zero = ...
    str2double(get(handles.edit_n_zero,'String'));

guidata(hObject, handles);

uiresume;


% --- Executes during object deletion, before destroying properties.
function datain_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to datain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function listbox_idata_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_idata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in listbox_idata.
function listbox_idata_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_idata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox_idata contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_idata

vars = evalin('base','who');
set(handles.listbox_idata,'String',vars)
guidata(hObject, handles);
% uiresume;


% --- Executes on button press in pushbutton_freq.
function pushbutton_freq_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    index_selected = get(handles.listbox_idata,'Value');
    var_list = get(handles.listbox_idata,'String');
    var_name = var_list{index_selected}; % Item selected in list box
    set(handles.text_freq,'String',var_name);
    guidata(hObject, handles)
%   uiresume

    
% --- Executes on button press in pushbutton_H_exp.
function pushbutton_H_exp_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_H_exp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    index_selected = get(handles.listbox_idata,'Value');
    var_list = get(handles.listbox_idata,'String');
    var_name = var_list{index_selected}; % Item selected in list box
    set(handles.text_H_exp,'String',var_name);
       
    [n_freq,n_io] = evalin('base',['size(',var_name,')']);
    set(handles.edit_n_o,'String',num2str(n_io));

	for ii=1:n_io,
		ind_def{ii} = num2str(0);
	end
	set(handles.edit_ind_SR_o,'String',ind_def);

	ind_def2{1} = num2str(0);
	set(handles.edit_ind_SR_i,'String',ind_def2);
        
    guidata(hObject, handles)
%   uiresume

    
% --- Executes on button press in pushbutton_Sig_H.
function pushbutton_Sig_H_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Sig_H (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    index_selected = get(handles.listbox_idata,'Value');
    var_list = get(handles.listbox_idata,'String');
    var_name = var_list{index_selected}; % Item selected in list box
    set(handles.text_Sig_H,'String',var_name);
    guidata(hObject, handles)
%   uiresume

    
% --- Executes during object creation, after setting all properties.
function edit_ind_SR_o_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ind_SR_o (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function edit_ind_SR_o_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ind_SR_o (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ind_SR_o as text
%        str2double(get(hObject,'String')) returns contents of edit_ind_SR_o as a double

guidata(hObject, handles);
% uiresume;


% --- Executes during object creation, after setting all properties.
function edit_ind_SR_i_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ind_SR_i (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function edit_ind_SR_i_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ind_SR_i (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ind_SR_i as text
%        str2double(get(hObject,'String')) returns contents of edit_ind_SR_i as a double

guidata(hObject, handles);
% uiresume;


% --- Executes during object creation, after setting all properties.
function edit_n_o_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_n_o (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function edit_n_o_Callback(hObject, eventdata, handles)
% hObject    handle to edit_n_o (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_n_o as text
%        str2double(get(hObject,'String')) returns contents of edit_n_o as a double

ind_no = str2double(get(hObject,'String'));
for ii=1:ind_no,
ind_def{ii} = num2str(0);
end
set(handles.edit_ind_SR_o,'String',ind_def);
guidata(hObject, handles)
% uiresume


% --- Executes during object creation, after setting all properties.
function edit_n_i_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_n_i (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function edit_n_i_Callback(hObject, eventdata, handles)
% hObject    handle to edit_n_i (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_n_i as text
%        str2double(get(hObject,'String')) returns contents of edit_n_i as a double

ind_ni = str2double(get(hObject,'String'));
for ii=1:ind_ni,
ind_def{ii} = num2str(0);
end
set(handles.edit_ind_SR_i,'String',ind_def);
guidata(hObject, handles);
% uiresume


% --- Executes during object creation, after setting all properties.
function edit_n_pole_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_n_pole (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function edit_n_pole_Callback(hObject, eventdata, handles)
% hObject    handle to edit_n_pole (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_n_pole as text
%        str2double(get(hObject,'String')) returns contents of edit_n_pole as a double

guidata(hObject, handles);
% uiresume;


% --- Executes during object creation, after setting all properties.
function edit_n_zero_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_n_zero (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function edit_n_zero_Callback(hObject, eventdata, handles)
% hObject    handle to edit_n_zero (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_n_zero as text
%        str2double(get(hObject,'String')) returns contents of edit_n_zero as a double

guidata(hObject, handles);
% uiresume;


% --- Executes on button press in pushbutton_coh.
function pushbutton_coh_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_coh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    index_selected = get(handles.listbox_idata,'Value');
    var_list = get(handles.listbox_idata,'String');
    var_name = var_list{index_selected}; % Item selected in list box
    set(handles.text_coh,'String',var_name);
    guidata(hObject, handles)
%   uiresume


% --- Executes on button press in pushbutton_Sig_H_inv.
function pushbutton_Sig_H_inv_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Sig_H_inv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    index_selected = get(handles.listbox_idata,'Value');
    var_list = get(handles.listbox_idata,'String');
    var_name = var_list{index_selected}; % Item selected in list box
    set(handles.text_Sig_H_inv,'String',var_name);
    guidata(hObject, handles)
%   uiresume


%
%	FINE
%
