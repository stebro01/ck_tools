function varargout = ck_rc2_multisubject_selection(varargin)
% CK_RC2_MULTISUBJECT_SELECTION M-file for ck_rc2_multisubject_selection.fig
%      CK_RC2_MULTISUBJECT_SELECTION, by itself, creates a new CK_RC2_MULTISUBJECT_SELECTION or raises the existing
%      singleton*.
%
%      H = CK_RC2_MULTISUBJECT_SELECTION returns the handle to a new CK_RC2_MULTISUBJECT_SELECTION or the handle to
%      the existing singleton*.
%
%      CK_RC2_MULTISUBJECT_SELECTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CK_RC2_MULTISUBJECT_SELECTION.M with the given input arguments.
%
%      CK_RC2_MULTISUBJECT_SELECTION('Property','Value',...) creates a new CK_RC2_MULTISUBJECT_SELECTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ck_rc2_multisubject_selection_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ck_rc2_multisubject_selection_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ck_rc2_multisubject_selection

% Last Modified by GUIDE v2.5 20-Jun-2017 22:46:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ck_rc2_multisubject_selection_OpeningFcn, ...
                   'gui_OutputFcn',  @ck_rc2_multisubject_selection_OutputFcn, ...
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


% --- Executes just before ck_rc2_multisubject_selection is made visible.
function ck_rc2_multisubject_selection_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ck_rc2_multisubject_selection (see VARARGIN)

% Choose default command line output for ck_rc2_multisubject_selection
handles.output = hObject;
handles.V = {};
set(handles.listbox1,'String',{});
set(handles.listbox1,'Value',1);
%set(handles.editsearchstring,'String','sw.*');
% Update handles structure
handles.startdir=pwd;
guidata(hObject, handles);

% UIWAIT makes ck_rc2_multisubject_selection wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ck_rc2_multisubject_selection_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

varargout{1} = handles.output1;
varargout{2} = handles.output2;
assignin('base','Vm',varargout{1});
assignin('base','Pm',varargout{2});
%varargout{1} = handles.output;
%assiginin('base','Pm',handles.output2);
guidata(hObject, handles);
close(handles.figure1);


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonadd.
function pushbuttonadd_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonadd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
list =get(handles.listbox1,'String');
assignin('base','addlist',list);
num_list = size(list,1)+1;
epifiles=spm_select([1 Inf],'image','Select EPI Files',[],...
    handles.startdir,get(handles.editsearchstring,'String'));

tmpdir=fileparts(epifiles(1,:));
tmploc=findstr(filesep,tmpdir);
handles.startdir=tmpdir(1:tmploc(1,size(tmploc,2))-1);
V = spm_vol(epifiles);

if size(handles.V,1)>0
    handles.V{size(handles.V,1)+1,1}=V;
    set(handles.listbox1,'Value',1);
else
    handles.V{1,1} = V;
end
% Subject-Name wird zum Speichern der Bilder benötigt
prompt = {'Enter Study/Proband name:'};
dlg_title = 'Input';
num_lines = 1;
def = {['Proband_',num2str(num_list)]};
answer = inputdlg(prompt,dlg_title,num_lines,def);
new_list = list;
new_list{num_list,1} = answer{1,1};
assignin('base','newlist',new_list);
set(handles.listbox1,'String',new_list);
%tmp = ck_get_nii_header_info(V(1,1).fname);
%voxel_dim = tmp.pixdim;
guidata(hObject,handles);

% --- Executes on button press in pushbuttonclearall.
function pushbuttonclearall_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonclearall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.listbox1,'String',{});
%set(handles.listbox1,'Value',0);
handles.V = {}; %rmfield(handles,'V');
guidata(hObject,handles);

% --- Executes on button press in pushbuttondeleteselected.
function pushbuttondeleteselected_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttondeleteselected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
idx = get(handles.listbox1,'Value');
%assignin('base','Vx',handles.V)

if idx>0
    list = get(handles.listbox1,'String');
    l = [1:idx-1,idx+1:size(list,1)];
        if size(l,2)==0
            new_list = {};
            rmfield(handles,'V');
        else
            new_list = list(l);
            handles.V=handles.V(l,1);
        end
    set(handles.listbox1,'String',new_list);
    if idx > length(new_list)
        set(handles.listbox1,'Value',length(new_list));
    end
end

guidata(hObject,handles);


% --- Executes on button press in pushbuttonok.
function pushbuttonok_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
numfiles=get(handles.editnumfiles,'String');
if strcmp('max',numfiles)==1
    handles.output1 = handles.V;
else
    n=str2num(numfiles);
    V=handles.V;
    for i=1:size(V,1)
        tmp=V{i,1}(1:n,1);
        V{i,1}=tmp;
    end
    handles.output1 = V;
end
handles.output2 = get(handles.listbox1,'String');
guidata(hObject,handles);
uiresume(handles.figure1);



function editsearchstring_Callback(hObject, eventdata, handles)
% hObject    handle to editsearchstring (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editsearchstring as text
%        str2double(get(hObject,'String')) returns contents of editsearchstring as a double


% --- Executes during object creation, after setting all properties.
function editsearchstring_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editsearchstring (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonsavelist.
function pushbuttonsavelist_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonsavelist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[f,p] = uiputfile('list.mat','Specify output file.');
Vm = handles.V;
Pm = get(handles.listbox1,'String');
save(fullfile(p,f),'Vm','Pm');

% --- Executes on button press in pushbuttonloadlist.
function pushbuttonloadlist_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonloadlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[f,p] = uigetfile();
load(fullfile(p,f));
handles.V = Vm;
set(handles.listbox1,'String',Pm);
guidata(hObject,handles);


% --- Executes on button press in pushbuttonmodifyallentries.
function pushbuttonmodifyallentries_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonmodifyallentries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Vm = handles.V;
Pm = get(handles.listbox1,'String');
assignin('base','Vm',Vm);
assignin('base','Pm',Pm);
mod_cell = inputdlg('Please enter modifier-string');
assignin('base','mod_cell',mod_cell);
mod = mod_cell{1,1};
for i=1:size(Vm,1)
    Pm{i,1} = [ mod Pm{i,1}];
    for j=1:size(Vm{i,1},1)
      [p,f,e] = fileparts(Vm{i,1}(j,1).fname);
      new_filename = fullfile(p,[mod f e]);
       Vm{i,1}(j,1).fname = new_filename;
    end
end
handles.V = Vm;
set(handles.listbox1,'String',Pm);
guidata(hObject,handles);



function editnumfiles_Callback(hObject, eventdata, handles)
% hObject    handle to editnumfiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editnumfiles as text
%        str2double(get(hObject,'String')) returns contents of editnumfiles as a double


% --- Executes during object creation, after setting all properties.
function editnumfiles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editnumfiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttondown.
function pushbuttondown_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttondown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
idx = get(handles.listbox1,'Value');
if idx<size(handles.V,1)
    TMP=handles.V(idx,1);
    handles.V(idx,1)=handles.V(idx+1,1);
    handles.V(idx+1,1)=TMP;
    Pm1 = get(handles.listbox1,'String');
    TMP=Pm1(idx,1);
    Pm1(idx,1)=Pm1(idx+1,1);
    Pm1(idx+1,1)=TMP;
    set(handles.listbox1,'String',[Pm1]);
    set(handles.listbox1,'Value',idx+1);
end
guidata(hObject,handles);

% --- Executes on button press in pushbuttonup.
function pushbuttonup_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
idx = get(handles.listbox1,'Value');
if idx>1
    TMP=handles.V(idx,1);
    handles.V(idx,1)=handles.V(idx-1,1);
    handles.V(idx-1,1)=TMP;
    Pm1 = get(handles.listbox1,'String');
    TMP=Pm1(idx,1);
    Pm1(idx,1)=Pm1(idx-1,1);
    Pm1(idx-1,1)=TMP;
    set(handles.listbox1,'String',[Pm1]);
    set(handles.listbox1,'Value',idx-1);
end
guidata(hObject,handles);
    

% --- Executes on button press in pushbuttonappendlist.
function pushbuttonappendlist_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonappendlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Vm1 = handles.V;
Pm1 = get(handles.listbox1,'String');
[f,p] = uigetfile();
load(fullfile(p,f));
handles.V = [Vm1;Vm];
set(handles.listbox1,'String',[Pm1;Pm]);
guidata(hObject,handles);


% --- Executes on button press in pushbuttonglr.
function pushbuttonglr_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonglr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% das ganze Verzeichnis wird durchsucht und die Subjects
% werden automatisch nach dem Template erstellt
% Namensgebung entspricht dem Verzeichnisnamen

path_org=pwd;
[p]=uigetdir();

%durchsuche rekursiv das Verzeichnis
%rrecursiv(p);

cd(p)
searchstring=get(handles.editsearchstring,'String');
befehl = ['dir /s /b ' searchstring];
fprintf ('>>Durchsuche Verzeichnis:...');
[status liste] = system(befehl);


% Übertrage die Namen in ein Cell-Array
o=1;
j=1; 
for i=1:size(liste,2) 
    if strcmp(liste(1,i),char(10)) 
        liste2{j,1}=liste(1,o:i-1); 
        o=i+1; 
        j=j+1; 
    end
end

fprintf('%g gefundene Dateien\n',length(liste2));

% Teile nach unterschiedlichen Verzeichnissen auf
[po f e] = fileparts(liste2{1,1});
j=1;
data=cell(1,1);
batch_idx=1;
hw = waitbar(0,['Durchsuche ' num2str(length(liste2)) ' Dateien']);

for i=1:size(liste2,1)
    
    if mod(i/length(liste2)*100,10) == 0
        waitbar(i/length(liste2),hw);
    end
    
    [p f e] = fileparts(liste2{i,1});
    if strcmp(p,po)
        % altes Verzeichnis
        data{j,1}=liste2{i,1};
        j=j+1;
    else
        % neues Verzeichnis
        data= get_sorted_files(data);
        assignin('base','data',data);
        
        % hier schlecht programmiert, das gleiche kommt unten nochmal
        % aber es funktioniert
list =get(handles.listbox1,'String');
assignin('base','addlist',list);
assignin('base','data',data);
num_list = size(list,1)+1;
%epifiles=spm_select([1 Inf],'image','Select EPI Files',[],...
%    handles.startdir,get(handles.editsearchstring,'String'));
epifiles=data;

tmpdir=fileparts(epifiles(1,:));
tmploc=findstr(filesep,tmpdir);
namedir=tmpdir(tmploc(1,size(tmploc,2))+1:end);
V = spm_vol(epifiles);

if size(handles.V,1)>0
    handles.V{size(handles.V,1)+1,1}=V;
    set(handles.listbox1,'Value',1);
else
    handles.V{1,1} = V;
end
% Subject-Name wird zum Speichern der Bilder benötigt
answer = namedir;
new_list = list;
new_list{num_list,1} = answer;
assignin('base','newlist',new_list);
set(handles.listbox1,'String',new_list);
        
        
        batch_idx=batch_idx+1;
        j=1;
        data=cell(1,1);
        po=p;
        data{j,1}=liste2{i,1};
        j=j+1;
    end
end
close(hw)

data= get_sorted_files(data);

list =get(handles.listbox1,'String');
assignin('base','addlist',list);
assignin('base','data',data);
num_list = size(list,1)+1;
%epifiles=spm_select([1 Inf],'image','Select EPI Files',[],...
%    handles.startdir,get(handles.editsearchstring,'String'));
epifiles=data;

tmpdir=fileparts(epifiles(1,:));
tmploc=findstr(filesep,tmpdir);
namedir=tmpdir(tmploc(1,size(tmploc,2))+1:end);
V = spm_vol(epifiles);

if size(handles.V,1)>0
    handles.V{size(handles.V,1)+1,1}=V;
    set(handles.listbox1,'Value',1);
else
    handles.V{1,1} = V;
end
% Subject-Name wird zum Speichern der Bilder benötigt
answer = namedir;
new_list = list;
new_list{num_list,1} = answer;
assignin('base','newlist',new_list);
set(handles.listbox1,'String',new_list);


cd(path_org)

guidata(hObject,handles);



function epifiles = get_sorted_files(data)


for i=1:size(data,1)
    sz(i)=size(data{i,1},2);
end

[ssz s_idx]=sort(sz);

data_org=data;
for i=1:size(data,1)
    data(i,1)=data_org(s_idx(i),1);
end
%epifiles=zeros(size(data,1),size(data{1,1},2)+2);
for i=1:size(data,1)
    epifiles(i,:)=[data{i,1} ',1'];
end
%epifiles=char(epifiles);
