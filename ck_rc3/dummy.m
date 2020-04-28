function varargout = ck_rc2(varargin)
% CK_RC2 MATLAB code for ck_rc2.fig
%      CK_RC2, by itself, creates a new CK_RC2 or raises the existing
%      singleton*.
%
%      H = CK_RC2 returns the handle to a new CK_RC2 or the handle to
%      the existing singleton*.
%
%      CK_RC2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CK_RC2.M with the given input arguments.
%
%      CK_RC2('Property','Value',...) creates a new CK_RC2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ck_rc2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ck_rc2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ck_rc2

% Last Modified by GUIDE v2.5 12-Mar-2018 11:59:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ck_rc2_OpeningFcn, ...
                   'gui_OutputFcn',  @ck_rc2_OutputFcn, ...
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


% --- Executes just before ck_rc2 is made visible.
function ck_rc2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ck_rc2 (see VARARGIN)

% Choose default command line output for ck_rc2
handles.output = hObject;
pushbuttonloadlast_Callback(hObject, eventdata, handles);
% axes(handles.axesstatus);
% xlim([0,1]);
% ylim([0,1]);
% xlabel = [];
% set(handles.axesstatus,'XTick',[]);
% set(handles.axesstatus,'YTick',[]);
%%%%%%%%%%%%%%%%%%%%
%construct waitbar
%grabs a handle 'f' to a figure. 
% f = gcf;  
% % Creates a WAITBAR in a new figure (by default) 
% h = waitbar(0,'Please wait...'); %   
% % The child of the waitbar is an axes object.  Grab the axes 
% % object 'c' and set its parent to be your figure f so that it now 
% % resides on figure f rather than on the old default figure. 
% c = get(h,'Children'); 
% set(c,'Parent',f);  % Set the position of the WAITBAR on your figure 
% set(c,'Units','Normalized','Position',[.05 .05 .8 .05]);  
% % Close the default figure 
% handles.waitbar = h;
% close(h);  

% The above steps only need to occur once to place the waitbar on the figure. 
% Now when you go to use the WAITBAR you can just call the waitbar with two 
% inputs: waitbar(x,h) where x is how much of the bar you want filled and f 
% is the handle to your figure. 




% Update handles structure
guidata(hObject, handles);
listboxupdate(hObject,handles);
%guidata(hObject, handles);

% UIWAIT makes ck_rc2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ck_rc2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function editnuminstances_Callback(hObject, eventdata, handles)
% hObject    handle to editnuminstances (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editnuminstances as text
%        str2double(get(hObject,'String')) returns contents of editnuminstances as a double


% --- Executes during object creation, after setting all properties.
function editnuminstances_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editnuminstances (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonnuminstancesless.
function pushbuttonnuminstancesless_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonnuminstancesless (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbuttonnuminstancesmore.
function pushbuttonnuminstancesmore_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonnuminstancesmore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbuttondataimport.
function pushbuttondataimport_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttondataimport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listboxupdate(hObject,handles);

choice = questdlg('What kind of data has to be analysed:', ...
    'choose kind of Data',...
    'MRI Data from 3D dicoms',...
    'MEG Data','Cancel','Cancel');
% Handle response
switch choice
    case 'MRI Data from 3D dicoms'
        if get(handles.checkboxgroup1,'Value')
            [Vmg1,Pmg1]=ck_rc2_multisubject_selection();
        end
        if get(handles.checkboxgroup2,'Value')
            [Vmg2,Pmg2]=ck_rc2_multisubject_selection();
        end
        if get(handles.checkboxgroup1,'Value')
            datadirg1 = get(handles.editdatadirg1,'String');
            ck_rc2_construct_data(datadirg1, Vmg1,Pmg1);
        end
        if get(handles.checkboxgroup2,'Value')
            datadirg2 = get(handles.editdatadirg2,'String');
            ck_rc2_construct_data(datadirg2, Vmg2,Pmg2);
        end
    case 'MRI Data from 4D dicoms'
        fprintf('noch nicht implementiert \n')
    case 'MEG Data'
        fprintf('noch nicht implementiert \n')    
    case 'Cancel'
end
parameterfile = fullfile(get(handles.editoutdir,'String'),'parameter.mat');
try
    load(parameterfile);
catch
    parameter = struct;
end

    TR = inputdlg('please enter TR','TR');
    parameter.TR = str2num(TR{1});
    lpf = inputdlg('please enter lp Filter Frequency','Frequency');
    parameter.lpf = str2num(lpf{1});
    hpf = inputdlg('please enter hp Filter Frequency','Frequency');
    parameter.hpf = str2num(hpf{1});
    parameter.bp_filter = [parameter.hpf parameter.lpf];
    save(parameterfile,'parameter');
    
listboxupdate(hObject,handles);




% --- Executes on selection change in listboxg1.
function listboxg1_Callback(hObject, eventdata, handles)
% hObject    handle to listboxg1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listboxg1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxg1


% --- Executes during object creation, after setting all properties.
function listboxg1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxg1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listboxg1 controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonparameter.
function pushbuttonparameter_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonparameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function editdatadirg1_Callback(hObject, eventdata, handles)
% hObject    handle to editdatadirg1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editdatadirg1 as text
%        str2double(get(hObject,'String')) returns contents of editdatadirg1 as a double


% --- Executes during object creation, after setting all properties.
function editdatadirg1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editdatadirg1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonupdate.
function pushbuttonupdate_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonupdate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listboxupdate(hObject,handles);


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbuttondefineclusters.
function pushbuttondefineclusters_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttondefineclusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% wechsel des Verzeichnesses da das Tool die Daten im gleichen Verzeichnis
% erwartet
cd(handles.datadirg1);
sb_roi_cluster_manager();
% falls im anderen Verzeichnis noch kein Cluster ist wird noch nach
% rueckfrage kopiert
% f1 = fullfile(get(handles.editdatadirg1,'String'),'Cluster.mat');
% f2 = fullfile(get(handles.editdatadirg2,'String'),'Cluster.mat');
% if exist(f2)~=2
%     copyfile(f1,f2);
% elseif exist(f2)==2
%     answer = questdlg('Soll Cluster.mat auch in den 2. datenordner kopiert werden?','overwrite?');
%     if strcmp(answer,'Yes')
%         copyfile(f1,f2);
%     end
% end
listboxupdate(hObject,handles);
cd(handles.outdir);



% --- Executes on button press in pushbuttondefinenetworks.
function pushbuttondefinenetworks_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttondefinenetworks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ck_rc2_cluster_grouping2(handles.outdir);

% --- Executes on button press in pushbuttonchangedirg1.
function pushbuttonchangedirg1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonchangedirg1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
p = uigetdir();
set(handles.editdatadirg1,'String',p);
guidata(hObject,handles);
listboxupdate(hObject,handles);


% --- Executes on button press in pushbuttonloadlast.
function pushbuttonloadlast_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonloadlast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[~,result] = dos('getmac');
mac = result(160:176);
try
    local_file = fullfile(fileparts(which('ck_rc2.m')),['lastsession_' mac '.mat']);
    load(local_file,'lastsession');
    set(handles.editoutdir,'String',lastsession.outdir);
    set(handles.editdatadirg1,'String',lastsession.datadirg1);
    set(handles.editdatadirg2,'String',lastsession.datadirg2);
    set(handles.editmultig1,'String',lastsession.multifileg1);
    set(handles.editmultig2,'String',lastsession.multifileg2);
    
    handles.outdir = lastsession.outdir;
    cd(handles.outdir);
catch
    fprintf('kein altes verzeichnis abgespeichert\n');
end
guidata(hObject,handles);
listboxupdate(hObject,handles);




% Update der listboxg1 mit den Variablen im aktuellen Verzeichnis
function listboxupdate(hObject,handles)


handles.outdir = get(handles.editoutdir,'String');
handles.datadirg1 = get(handles.editdatadirg1,'String');
handles.datadirg2 = get(handles.editdatadirg2,'String');
handles.multifileg1 = get(handles.editmultig1,'String');
handles.multifileg2 = get(handles.editmultig2,'String');

%set(handles.listboxg1,'Value',1);
%set(handles.listboxg2,'Value',1);

g1 = dir([handles.datadirg1 filesep '*.mat']);
g2 = dir([handles.datadirg2 filesep '*.mat']);
if get(handles.checkboxfilter,'Value')
    filterstr = get(handles.editfilter,'String');
    g1 = get_filtered_list(g1,filterstr);
    g2 = get_filtered_list(g2,filterstr);
end

assignin('base','g1',g1);

listg1{1}='empty';
idxg1 = 1;
listg2{1}='empty';
idxg2 = 1;
%Val{1} ='empty';
preselect = get(handles.checkboxpreselect,'Value');
preselectstr = get(handles.editpreselect,'String');
Valg1(1)=1;
Valg2(1)=1;
for i=1:size(g1,1)
    listg1{i} = g1(i,1).name;
    if preselect
        if length(preselectstr)<=length(g1(i,1).name)
           
        if strcmp(g1(i,1).name(1:length(preselectstr)),preselectstr)
            Valg1(idxg1) = i;
            idxg1 = idxg1 + 1;
        end
        end
    end
end
for i=1:size(g2,1)
    listg2{i} = g2(i,1).name;
    if preselect
        if length(preselectstr)<=length(g2(i,1).name)
           
        if strcmp(g2(i,1).name(1:length(preselectstr)),preselectstr)
                    Valg2(idxg2) = i;
            idxg2 = idxg2 + 1;
        end
        end
    end
end
assignin('base','listg1',listg1);
set(handles.listboxg1,'String',listg1);
set(handles.listboxg2,'String',listg2);
% tmp=get(handles.listboxg1,'Value');
% assignin('base','tmp',tmp);
if idxg1>1
    set(handles.listboxg1,'Value',Valg1);
else 
    set(handles.listboxg1,'Value',1);
end
if idxg2>1
    set(handles.listboxg2,'Value',Valg2);
else
    set(handles.listboxg2,'Value',1);
end

%save(fullfile(fileparts(which('ck_rc2.m')),'lastprocessingdir.mat'),'lastprocessingdir');
% speichere REchnerspezifisch unter Hilfe der MAC Adresse des REchners
[~,result] = dos('getmac');
mac = result(160:176);
try
    load(fileparts(which(which('ck_rc2.m')),['lastsession_' mac '.mat']),'lastsession');

    cd(handles.outdir);
catch
    lastsession.outdir = handles.outdir;
    cd(handles.outdir);
end
    lastsession.outdir = handles.outdir;
    lastsession.datadirg1 = handles.datadirg1;
    lastsession.datadirg2 = handles.datadirg2;
    lastsession.multifileg1 = handles.multifileg1;
    lastsession.multifileg2 = handles.multifileg2;
savefile = fullfile(fileparts(which('ck_rc2.m')),['lastsession_' mac '.mat']);
save(savefile,'lastsession');
parameterfile = fullfile(handles.outdir,'parameter.mat');
try  
    load(parameterfile);
catch
    parameter = struct;
end
parameter.outdir = handles.outdir;
parameter.datadirg1 = handles.datadirg1;
parameter.datadirg2 = handles.datadirg2;
try
save(parameterfile,'parameter');
catch
end
guidata(hObject,handles);


function g = get_filtered_list(g1,filterstr)
% sortiert alle eintraege aus die nicht dem filterstr entsprechen
idx = 0;
for i=1:length(g1)
    str = g1(i).name;
    k = strfind(str,filterstr);
    if length(k)>0
        idx = idx + 1;
        g(idx,1)=g1(i);
    end
end

    




% --- Executes on button press in pushbuttondataimportinfo.
function pushbuttondataimportinfo_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttondataimportinfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbuttonfillclusterwithdata.
function pushbuttonfillclusterwithdata_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonfillclusterwithdata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%listboxupdate(hObject,handles);

reindexing = handles.checkboxreindexing.Value;
[filenamesg1, filenamesg2] = get_filenames(handles,'Data_');
%[G1, G2] = get_all_selected_filenames(handles);

if  get(handles.checkboxfillclusterwithdatag1,'Value')
    ck_rc2_fill_all_clusters(get(handles.editdatadirg1,'String'),filenamesg1,'Cluster',reindexing);
    if handles.checkboxincludeAAL.Value
        ck_rc2_fill_all_clusters(get(handles.editdatadirg1,'String'),filenamesg1,'AALClustAAL',reindexing);
    end
    
end
if  get(handles.checkboxfillclusterwithdatag2,'Value')
    ck_rc2_fill_all_clusters(get(handles.editdatadirg2,'String'),filenamesg2,'Cluster',reindexing);
    if handles.checkboxincludeAAL.Value
        ck_rc2_fill_all_clusters(get(handles.editdatadirg2,'String'),filenamesg2,'AALClustAAL',reindexing);
    end
end
listboxupdate(hObject,handles);
% 

% --- Executes on button press in pushbuttoninfoparameter.
function pushbuttoninfoparameter_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttoninfoparameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbuttoninfofillcluster.
function pushbuttoninfofillcluster_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttoninfofillcluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

web  ck_rc2_info_fillcluster.htm -browser 

% --- Executes on button press in pushbuttoninfoestimateconnectivity.
function pushbuttoninfoestimateconnectivity_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttoninfoestimateconnectivity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton19.
function pushbutton19_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbuttondefinenetworkmulti.
function pushbuttondefinenetworkmulti_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttondefinenetworkmulti (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ck_rc2_cluster_grouping_multi( get(handles.editoutdir,'String') )

% 
% ck_rc2_cluster_grouping_multi( get(handles.editdatadirg1,'String') )
% f1 = fullfile(get(handles.editdatadirg1,'String'),'Cluster.mat');
% f2 = fullfile(get(handles.editdatadirg2,'String'),'Cluster.mat');
% if exist(f2)~=2
%     copyfile(f1,f2);
% elseif exist(f2)==2
%     answer = questdlg('Soll Cluster.mat auch in den 2. datenordner kopiert werden?','overwrite?');
%     if strcmp(answer,'Yes')
%         copyfile(f1,f2);
%     end
% end
% 
% answer = questdlg('Soll die Netzwerkarchitektur auch in alle anderen Cluster uebertragen werden?','overwrite?');
% if strcmp(answer,'Yes')
%     update_network(get(handles.editdatadirg1,'String'));
%     update_network(get(handles.editdatadirg2,'String'));
%     fprintf('all Cluster files updated by new networks');
% end

function update_network(dirname)

%lade Clusterstruktur mit den aktualisieren Netzwerkeintraegen
load(fullfile(dirname,'Cluster.mat'));
Cluster_org = Cluster;
t = dir(fullfile(dirname,'*Cluster*.mat'));

for i=1:length(t)
    load(fullfile(dirname,t(i).name));
    for j=1:length(Cluster)
        Cluster{j}.net = Cluster_org{i}.net;
        if isfield(Cluster{j},'num_global_networks')
            Cluster{j}.num_global_networks = Cluster_org{i}.num_global_networks;
        end
        if isfield(Cluster{j},'num_local_networks')
            Cluster{j}.num_local_networks = Cluster_org{i}.num_local_networks;
        end
    end
    save(fullfile(dirname,t(i).name),'Cluster','-v7.3');
end


function editdatadirg2_Callback(hObject, eventdata, handles)
% hObject    handle to editdatadirg2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editdatadirg2 as text
%        str2double(get(hObject,'String')) returns contents of editdatadirg2 as a double


% --- Executes during object creation, after setting all properties.
function editdatadirg2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editdatadirg2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobuttong1.
function radiobuttong1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttong1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttong1


% --- Executes on selection change in listboxg2.
function listboxg2_Callback(hObject, eventdata, handles)
% hObject    handle to listboxg2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listboxg2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxg2


% --- Executes during object creation, after setting all properties.
function listboxg2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxg2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobuttong2.
function radiobuttong2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttong2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttong2


% --- Executes on button press in pushbuttonload.
function pushbuttonload_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonload (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
p = uigetdir();
parameterfile = fullfile(p,'parameter.mat');
if ~exist(parameterfile)
    warndlg('no parameter.mat file in choosen directory');
    return;
end
load(parameterfile);
parameter.datadirg1
set(handles.editoutdir,'String',parameter.outdir);
set(handles.editdatadirg1,'String',parameter.datadirg1);
set(handles.editdatadirg2,'String',parameter.datadirg2);

guidata(hObject,handles);
listboxupdate(hObject,handles);



function editfilter_Callback(hObject, eventdata, handles)
% hObject    handle to editfilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editfilter as text
%        str2double(get(hObject,'String')) returns contents of editfilter as a double


% --- Executes during object creation, after setting all properties.
function editfilter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editfilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxfilter.
function checkboxfilter_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxfilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxfilter


% --- Executes on button press in pushbuttonchangedirg2.
function pushbuttonchangedirg2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonchangedirg2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
p = uigetdir();
set(handles.editdatadirg2,'String',p);
guidata(hObject,handles);
listboxupdate(hObject,handles);



function editoutdir_Callback(hObject, eventdata, handles)
% hObject    handle to editoutdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editoutdir as text
%        str2double(get(hObject,'String')) returns contents of editoutdir as a double


% --- Executes during object creation, after setting all properties.
function editoutdir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editoutdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonchangediroutdir.
function pushbuttonchangediroutdir_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonchangediroutdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
p = uigetdir();
cd(p);
handles.outdir = p;
set(handles.editoutdir,'String',p);
guidata(hObject,handles);
listboxupdate(hObject,handles);


% --- Executes on button press in radiobuttontakeall.
function radiobuttontakeall_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttontakeall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttontakeall
if get(hObject,'Value')
    set(handles.radiobuttontakeselected,'Value',0);
else
    set(handles.radiobuttontakeselected,'Value',1);
end

% --- Executes on button press in radiobuttontakeselected.
function radiobuttontakeselected_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttontakeselected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttontakeselected
if get(hObject,'Value')
    set(handles.radiobuttontakeall,'Value',0);
else
     set(handles.radiobuttontakeall,'Value',1);
end

% --- Executes on button press in checkboxgroup1.
function checkboxgroup1_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxgroup1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxgroup1


% --- Executes on button press in checkboxgroup2.
function checkboxgroup2_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxgroup2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxgroup2


% --- Executes on button press in pushbuttonsave.
function pushbuttonsave_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonsave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listboxupdate(hObject,handles);



function editpreselect_Callback(hObject, eventdata, handles)
% hObject    handle to editpreselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editpreselect as text
%        str2double(get(hObject,'String')) returns contents of editpreselect as a double


% --- Executes during object creation, after setting all properties.
function editpreselect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editpreselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxpreselect.
function checkboxpreselect_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxpreselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxpreselect


% --- Executes on button press in pushbuttong2tog1.
function pushbuttong2tog1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttong2tog1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listboxupdate(hObject,handles);
idx = get(handles.listboxg2,'Value');
list = get(handles.listboxg2,'String');
filenames = list(idx);
assignin('base','filenames',filenames);
for i=1:length(filenames)
    f1 = fullfile(get(handles.editdatadirg2,'String'),filenames{i});
    f2 = fullfile(get(handles.editdatadirg1,'String'),filenames{i});
    copyfile(f1,f2);
end
listboxupdate(hObject,handles);

% --- Executes on button press in pushbuttong1tog2.
function pushbuttong1tog2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttong1tog2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listboxupdate(hObject,handles);
idx = get(handles.listboxg1,'Value');
list = get(handles.listboxg1,'String');
filenames = list(idx);
assignin('base','filenames',filenames);
for i=1:length(filenames)
    f1 = fullfile(get(handles.editdatadirg1,'String'),filenames{i});
    f2 = fullfile(get(handles.editdatadirg2,'String'),filenames{i});
    copyfile(f1,f2);
end
listboxupdate(hObject,handles);
        


% --- Executes on button press in pushbuttoninspectdata.
function pushbuttoninspectdata_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttoninspectdata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Modified CK 28.11.2017
% [filenamesg1, filenamesg2] = get_filenames(handles,'Cluster_');
% d1 = get(handles.editdatadirg1,'String');
% d2 = get(handles.editdatadirg2,'String');
% for i=1:length(filenamesg1)
%     ffilenamesg1{i}=fullfile(d1,filenamesg1{i});
% end
% for i=1:length(filenamesg2)
%     ffilenamesg2{i}=fullfile(d2,filenamesg2{i});
% end
[ffilenamesg1,ffilenamesg2] = get_all_selected_filenames(handles,'part');

% end of modification CK
load(fullfile(handles.outdir,'parameter.mat'));
parameter.tmpffilenamesg1 = ffilenamesg1;
parameter.tmpffilenamesg2 = ffilenamesg2;
save(fullfile(handles.outdir,'parameter.mat'),'parameter');
%ck_rc2_inspect_data(handles.outdir,ffilenamesg1,ffilenamesg2);
tmp = pwd;
cd(handles.outdir);
ck_rc2_inspect_data_control()
cd(tmp);


% --- Executes on button press in checkboxincludeAAL.
function checkboxincludeAAL_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxincludeAAL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxincludeAAL
% copiere den AALCluster in beide Verzeichnisse falls noch nicht vorhanden
%load('AALClusterAAL.mat')


% --- Executes on button press in pushbuttonreadmeincludeAAL.
function pushbuttonreadmeincludeAAL_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonreadmeincludeAAL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
readme{1,1} = 'Include AAL\n';
readme{2,1} = '-----------------------------------------------------------\n';
readme{3,1} = 'Es soll zusaetzlich eine Berechnung fuer eine Clusterstruktur\n';
readme{4,1} = 'erfolgen die dem AAL Atlas entspricht\n';
readme{5,1} = 'dazu muss in beiden Datenverzeichnissen das matfile "AALClustAAL" \n';
readme{6,1} = 'vorhanden sein. Ein Template hierfuer befindet sich im ck_rc2 Ordner \n';
readme{7,1} = 'das muss dann aber noch von der Aufloesung her angepasst werden \n';
readme{8,1} = '(adapt resolution in Define Clusters \n';
readme{9,1} = 'Nur wenn dieses File vorhanden ist kann spaeter die Optionen \n';
readme{10,1} = 'AAL2AAL und Cluster2AAL in der connectivite estimation genutzt werden \n';

txt = '';
for i=1:size(readme,1)
    txt = [txt readme{i,1}]; %#ok<AGROW>
end
helpdlg(sprintf(txt));


% --- Executes on button press in pushbuttondeleteg1.
function pushbuttondeleteg1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttondeleteg1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
answer = questdlg('do you really want to delete these files?');
if strcmp(answer,'Yes')

d1 = [get(handles.editdatadirg1,'String'),filesep'];
listg1 = handles.listboxg1.String(handles.listboxg1.Value);
listg1 = strcat(d1,listg1);
for i=1:length(listg1)
    delete(listg1{i});
end
end
listboxupdate(hObject,handles);





% --- Executes on button press in pushbuttondeleteg2.
function pushbuttondeleteg2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttondeleteg2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
answer = questdlg('do you really want to delete these files?');
if strcmp(answer,'Yes')
d2 = [get(handles.editdatadirg2,'String'),filesep'];
listg2 = handles.listboxg2.String(handles.listboxg2.Value);
listg2 = strcat(d2,listg2);
for i=1:length(listg2)
    delete(listg2{i});
end
end
listboxupdate(hObject,handles);


% --- Executes on button press in checkboxreindexing.
function checkboxreindexing_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxreindexing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxreindexing


% --- Executes on button press in pushbuttonreducedata.
function pushbuttonreducedata_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonreducedata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


[listg1,listg2] = get_all_selected_filenames(handles,'full');
prefix = handles.editreduceclusternewname.String;
reduce_Cluster(handles,listg1,prefix);
reduce_Cluster(handles,listg2,prefix);


%guidata(hObject,handles);
listboxupdate(hObject,handles);


function reduce_Cluster(handles,filelist,prefix)
% entfernt die Y Daten und speichert die Cluster unter neuem Namen ab

for i=1:length(filelist)
    fprintf('reduce_Cluster %.0d / %.0d subjects\n',i,length(filelist));
    file = filelist{i};
    [p,f,e] = fileparts(file);
    newfilename = fullfile(p,[prefix, f, e]);
    load(file); % load Cluster
    for j=1:length(Cluster)
        
        if isfield(Cluster{j},'Y') && handles.checkboxredY.Value
            Cluster{j} = rmfield(Cluster{j},'Y');
        end
        if isfield(Cluster{j},'XYZ') && handles.checkboxredXYZ.Value
        Cluster{j} = rmfield(Cluster{j},'XYZ');
        end
        if isfield(Cluster{j},'Z') && handles.checkboxredZ.Value
            Cluster{j} = rmfield(Cluster{j},'Z');
        end 
           
        if isfield(Cluster{j},'XYZ_ind') && handles.checkboxredXYZ_ind.Value
            Cluster{j} = rmfield(Cluster{j},'XYZ_ind');
        end
        if isfield(Cluster{j},'XZY_ind_all')&& handles.checkboxredXYZ_ind_all.Value
            Cluster{j} = rmfield(Cluster{j},'XZY_ind_all');
        end
    end
    save(newfilename,'Cluster');
    
end





function [listg1,listg2] = get_all_selected_filenames(handles,string)
% gibt alle selektierten filenamen zuruck
% string kennzeichnet die Option
% string = 'full' ... voller filename mit path
% string = 'part' ... partieller filename ohne path

d1 = [get(handles.editdatadirg1,'String'),filesep'];
d2 = [get(handles.editdatadirg2,'String'),filesep'];

listg1 = handles.listboxg1.String(handles.listboxg1.Value);
listg2 = handles.listboxg2.String(handles.listboxg2.Value);

if strcmp(string,'full')
    listg1 = strcat(d1,listg1);
    listg2 = strcat(d2,listg2);
end


function [filenamesg1,filenamesg2] = get_filenames(handles,pattern)

filenamesg1={};
filenamesg2={};
% 
% if get(handles.radiobuttontakeall,'Value') && get(handles.checkboxgroup1,'Value')
%     % alle Datafiles
%     idx = 0;
%     listg1tmp=get(handles.listboxg1,'String');
%     for i=1:size(listg1tmp,1)
%         if strcmp(listg1tmp{i}(1:length(pattern)),pattern)
%             idx = idx + 1;
%             filenamesg1{idx}=listg1tmp{i};
%         end
%     end
%         
% else
if get(handles.checkboxgroup1,'Value')
    idxg1 = get(handles.listboxg1,'Value');
    listg1 = get(handles.listboxg1,'String');
    filenamesg1 = listg1(idxg1);
 end


% erneut fuer datadirectory 2
% if get(handles.radiobuttontakeall,'Value') && get(handles.checkboxgroup2,'Value')
%     % alle Datafiles
%     idx = 0;
%     listg2tmp=get(handles.listboxg2,'String');
%     for i=1:size(listg2tmp,1)
%         if strcmp(listg2tmp{i}(1:length(pattern)),pattern)
%             idx = idx + 1;
%             filenamesg2{idx}=listg2tmp{i};
%         end
%     end
%         
% else
if get(handles.checkboxgroup2,'Value')
    idxg2 = get(handles.listboxg2,'Value');
    listg2 = get(handles.listboxg2,'String');
    filenamesg2 = listg2(idxg2);
 end


% --- Executes on button press in pushbuttondeselectg1.
function pushbuttondeselectg1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttondeselectg1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.listboxg1,'value',[]); 
guidata(hObject,handles);

% --- Executes on button press in pushbuttondeselectg2.
function pushbuttondeselectg2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttondeselectg2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.listboxg2,'value',[]); 
guidata(hObject,handles);


function editreduceclusternewname_Callback(hObject, eventdata, handles)
% hObject    handle to editreduceclusternewname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editreduceclusternewname as text
%        str2double(get(hObject,'String')) returns contents of editreduceclusternewname as a double


% --- Executes during object creation, after setting all properties.
function editreduceclusternewname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editreduceclusternewname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxfillclusterwithdatag1.
function checkboxfillclusterwithdatag1_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxfillclusterwithdatag1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxfillclusterwithdatag1


% --- Executes on button press in checkboxfillclusterwithdatag2.
function checkboxfillclusterwithdatag2_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxfillclusterwithdatag2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxfillclusterwithdatag2


% --- Executes on button press in checkboxestimateconnectivityg1.
function checkboxestimateconnectivityg1_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxestimateconnectivityg1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxestimateconnectivityg1


% --- Executes on button press in checkboxestimateconnectivityg2.
function checkboxestimateconnectivityg2_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxestimateconnectivityg2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxestimateconnectivityg2


% --- Executes on button press in checkboxredY.
function checkboxredY_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxredY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxredY


% --- Executes on button press in checkboxredXYZ.
function checkboxredXYZ_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxredXYZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxredXYZ


% --- Executes on button press in checkboxredZ.
function checkboxredZ_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxredZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxredZ


% --- Executes on button press in checkboxredXYZ_ind.
function checkboxredXYZ_ind_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxredXYZ_ind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxredXYZ_ind


% --- Executes on button press in checkboxredXZY_ind_all.
function checkboxredXZY_ind_all_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxredXZY_ind_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxredXZY_ind_all


% --- Executes on button press in pushbuttonreadmereducedata.
function pushbuttonreadmereducedata_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonreadmereducedata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
readme{1,1} = 'reduce Data\n';
readme{2,1} = '-----------------------------------------------------------\n';
readme{3,1} = 'Wenn sehr grosse Clusterstrukturen insbesondere der AAL Atlas\n';
readme{4,1} = 'benutzt werden, werden die Clusterfiles der Einzelprobanden\n';
readme{5,1} = 'so gross das das einladen der Clusterstrukturen viel Zeit und  \n';
readme{6,1} = 'Speicher frist.\n';
readme{7,1} = 'Fuer das Berechnen und visulisieren der Konnektivitaeten werden \n';
readme{8,1} = 'jedoch zumeist die Ursprungsdaten nicht mehr benoetigt sondern\n';
readme{9,1} = 'mehr z.B. Y_mean \n';
readme{10,1} = 'Mit reduce Data koennen daher bestimmte Felder in der Clusterstruktur \n';
readme{11,1} = 'geloescht werden und der Cluster mit einem prefix neu abgespeichert werden \n';

txt = '';
for i=1:size(readme,1)
    txt = [txt readme{i,1}]; %#ok<AGROW>
end
helpdlg(sprintf(txt));


% --- Executes on button press in pushbuttontest.
function pushbuttontest_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttontest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% [g1,g2] = get_all_selected_filenames(handles,'full');
% for i=1:length(g2)
%     load(g2{i});
%     set_status(handles,'working',i/length(g2));
% end
steps = 10; 
rectangle('Position',[0,0,1,1],'FaceColor','white')
for step = 1:steps     
  % computations take place here     
%  waitbar(step / steps,handles.waitbar)     
  set_status('test',step / steps)     
  pause(.34) 
end 
rectangle('Position',[0,0,1,1],'FaceColor','white')



function set_status(description, p)

%axes(handles.axesstatus);
rectangle('Position',[0,0,p,1],'FaceColor','r')
title(description)
%daspect([1,1,1])


% --- Executes on button press in pushbuttoncspd.
function pushbuttoncspd_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttoncspd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ck_rc2_getting_and_structuring_preprocessed_data();


% --- Executes on button press in pushbuttondeletenetworkdef.
function pushbuttondeletenetworkdef_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttondeletenetworkdef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h=waitbar(0,'Please wait ...');
try
    f = fullfile(handles.editoutdir.String,'Cluster.mat');
load(f);
Cluster = del_network(Cluster);
    save(f,'Cluster');
 catch
     fprintf('catch during puschbuttondeletnetworkdef in erster \n');
 end
waitbar(1 / 3);
try
load(fullfile(handles.editdatadirg1.String,'Cluster.mat'));
Cluster = del_network(Cluster);
save(fullfile(handles.editdatadirg1.String,'Cluster.mat'),'Cluster');
waitbar(2/3);
load(fullfile(handles.editdatadirg2.String,'Cluster.mat'));
Cluster = del_network(Cluster);
save(fullfile(handles.editdatadirg2.String,'Cluster.mat'),'Cluster');
catch
    fprintf('catch during puschbuttondeletnetworkdef in zweiter \n');
end
waitbar(3/3)
close(h);

function C = del_network(C)
for i=1:length(C)
    if isfield(C{i},'net')
        C{i} = rmfield(C{i},'net');
    end
    if isfield(C{i},'group')
        C{i} = rmfield(C{i},'group');
    end
    if isfield(C{i},'groupname')
        C{i} = rmfield(C{i},'groupname');
    end
    if isfield(C{i},'groupm')
        C{i} = rmfield(C{i},'groupm');
    end
    if isfield(C{i},'groupnamem')
        C{i} = rmfield(C{i},'groupnamem');
    end
    if isfield(C{i},'num_global_networks')
        C{i} = rmfield(C{i},'num_global_networks');
    end
    if isfield(C{i},'num_local_networks')
        C{i} = rmfield(C{i},'num_local_networks');
    end
end


% --- Executes on button press in pushbuttonfillclusterfromniftifiles.
function pushbuttonfillclusterfromniftifiles_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonfillclusterfromniftifiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
m1file = get(handles.editmultig1,'String');
load(m1file);
Vmg1 = Vm;
Pmg1 = Pm;
m2file = get(handles.editmultig2,'String');
load(m2file);
Vmg2 = Vm;
Pmg2 = Pm;
%options = ck_rc2_fill_all_clusters_from_niftifiles_options(Vmg1,Pmg1,Vmg2,Pmg2);
%options = ck_rc2_fill_all_clusters_from_niftifiles_options();
parameterfile = fullfile(get(handles.editoutdir,'String'),'parameter.mat');
hApp = ck_rc2_fill_all_clusters_from_niftifiles_options;
options = MySetupFunction(hApp, Vmg1,Pmg1,Vmg2,Pmg2,parameterfile);
%[Vmg1,Pmg1]=ck_rc2_multisubject_selection();
%[Vmg2,Pmg2]=ck_rc2_multisubject_selection();
datadirg1 = get(handles.editdatadirg1,'String');
datadirg2 = get(handles.editdatadirg2,'String');
clusterdir = get(handles.editoutdir,'String');
if handles.checkboxfillclusterwithdatag1.Value
    ck_rc2_fill_all_clusters_from_niftifiles(clusterdir,datadirg1,Vmg1,Pmg1,options);
end
if handles.checkboxfillclusterwithdatag2.Value
    ck_rc2_fill_all_clusters_from_niftifiles(clusterdir,datadirg2,Vmg2,Pmg2,options);
end


function editmultig1_Callback(hObject, eventdata, handles)
% hObject    handle to editmultig1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editmultig1 as text
%        str2double(get(hObject,'String')) returns contents of editmultig1 as a double


% --- Executes during object creation, after setting all properties.
function editmultig1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editmultig1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonselectmultig1.
function pushbuttonselectmultig1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonselectmultig1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[f,p] = uigetfile();
set(handles.editmultig1,'String',fullfile(p,f));



function editmultig2_Callback(hObject, eventdata, handles)
% hObject    handle to editmultig2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editmultig2 as text
%        str2double(get(hObject,'String')) returns contents of editmultig2 as a double


% --- Executes during object creation, after setting all properties.
function editmultig2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editmultig2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonselectmultig2.
function pushbuttonselectmultig2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonselectmultig2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[f,p] = uigetfile();
set(handles.editmultig2,'String',fullfile(p,f));

% --- Executes on button press in pushbuttoncreatemulti.
function pushbuttoncreatemulti_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttoncreatemulti (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[V,P] = ck_rc2_multisubject_selection();



% --- Executes on button press in pushbuttonestimateconnectivity.
function pushbuttonestimateconnectivity_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonestimateconnectivity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%filenames = cell{length(idx),1};


%listboxupdate(hObject,handles);

%[filenamesg1,filenamesg2] = get_all_selected_fullfilenames(handles);
[filenamesg1,filenamesg2] = get_all_selected_filenames(handles,'part');
assignin('base','file',filenamesg1);
%[filenamesg1, filenamesg2] = get_filenames(handles,'Cluster_');
outdir = get(handles.editoutdir,'String');
load(fullfile(outdir,'parameter.mat'));
options = ck_rc2_estimate_connectivity_options();
options.parameter = parameter;
%if  get(handles.checkboxestimateconnectivityg1,'Value')

if  length(filenamesg1)>0
   ck_rc2_estimate_connectivity(handles.outdir,handles.datadirg1,filenamesg1,options);
end
%if  get(handles.checkboxestimateconnectivityg2,'Value')
if  length(filenamesg2)>0
   ck_rc2_estimate_connectivity(handles.outdir,handles.datadirg2,filenamesg2,options);
end



listboxupdate(hObject,handles);




% --- Executes on button press in pushbuttonvisualiseandtest.
function pushbuttonvisualiseandtest_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonvisualiseandtest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

assignin('base','handles',handles);
handles.outdir = get(handles.editoutdir,'String');
handles.datadirg1 = get(handles.editdatadirg1,'String');
handles.datadirg2 = get(handles.editdatadirg2,'String');

str={'one group Cluster2Cluster';'two groups Cluster2Cluster';'AAL2AAL';...
    'Voxel_in_AAL2AAL';'ICA';'one group mtrail granger C2C';'one group granger C2C'...
    ;'three groups C2C';'Cluster2Cluster side differences';'one group behavioral'};
[s,v] = listdlg('PromptString','Please choose:',...
                'SelectionMode','single',...
                'ListString',str);
if v==0
    return
end
% Handle response
switch s
    case 1
        ck_rc_visualise_C2C_one_group_control(handles.outdir);
        warndlg('currently not adaptaed from ck_resting_conn()');
%        break   
    case 2
        ck_rc2_visualise_C2C_control2({handles.outdir, handles.datadirg1, handles.datadirg2});
%        break   
    case 3
        ck_rc_visualise_AAL2AAL_control(handles.outdir);
        warndlg('currently not adaptaed from ck_resting_conn()');
%        break
    case 4
        ck_rc_visualise_VoxelinAAL2AAL_control(handles.outdir);
        warndlg('currently not adaptaed from ck_resting_conn()');
    case 5
        ck_rc_visualise_ICA_control(handles.outdir);
        warndlg('currently not adaptaed from ck_resting_conn()');
    case 6
        ck_rc_visualize_C2C_granger_one_group_mtrail_control(handles.outdir);
        warndlg('currently not adaptaed from ck_resting_conn()');
    case 7
        ck_rc_visualize_C2C_granger_one_group_control(handles.outdir);
        warndlg('currently not adaptaed from ck_resting_conn()');
    case 8
        ck_rc_visualize_C2C_three_groups_control(handles.outdir);
        warndlg('currently not adaptaed from ck_resting_conn()');
    case 9
        ck_rc_visualize_C2C_side_differences_control(handles.outdir);
        warndlg('currently not adaptaed from ck_resting_conn()');
    case 10
        ck_rc2_visualize_C2behavior_control(handles.outdir);
    otherwise
        % do nothing
end

