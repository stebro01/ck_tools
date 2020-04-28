function [ output_args ] = ck_rc3_define_network( outdir )
%UNTITLED GUI for selection of groups of Clusters
%   Detailed explanation goes here

% Positioniere das neue Fenster wenn es passt neben das bestehende
% sonst in die obere linke Ecke
%pz = get(pfh,'Position');
% GUI datan
handles = struct;
handles.output = struct;
handles.output.cancel = 0;
handles.outdir = outdir;

try
    load(fullfile(handles.outdir,'Cluster.mat'));
catch 
    warndlg('no Cluster structure available in current directory','Warning');
    handles.output.cancel = 1;
    guidata(hObject,handles);
    return;
end

num_cluster = length(Cluster);

sz = get(0,'ScreenSize');
width = 1500;
maxheight = sz(4)-200;
height = 100 + num_cluster *30 + 40;
stop = 0;
if height> (sz(4)-100)
    height = maxheight;
    fprintf('Clusteranzahln nicht mehr auf einer Seite darstellbar\n');
    for i=1:num_cluster
        h = 100 + i *30 + 40;
        if h>maxheight && stop==0
            num_vis = i-1;
            stop = 1;
        end
    end
else
    num_vis = num_cluster;
end
f = figure('Visible','off','Position',[50,50,width,height],...
    'MenuBar','none',...
    'Resize','off',...
    'Toolbar','none');



%Define Buttons
bw = 90;
bh = 20;
%Define Editfields
ew = 80;
eh = 20;




% converting Cluster
for i=1:length(Cluster)
    Clusternames{i,1} = Cluster{i}.name;
    middle = floor(length(Cluster{i}.XYZ)/2);
    try
    Clustercoor{i,1} = num2str(round(Cluster{i}.XYZ(:,middle)'));
    catch ME
        fprintf('error in ck_rc2_cluster_grouping_multi\n');
        fprintf('"Clustercoor{i,1} = num2str(round(Cluster{i}.XYZ(:,middle)));" \n');
        fprintf('i = %.0f   middle = %.0f\n');
        assignin('base','Cluster',Cluster);
        rethrow ME
    end
      
    if isfield(Cluster{i},'group')
        Clustergroup{i,1} = Cluster{i}.group;
    else
        Clustergroup{i,1} = 1;
    end
    if isfield(Cluster{i},'groupname')
        Clustergroupname{i,1} = Cluster{i}.groupname;
    else
        Clustergroupname{i,1} = 'noname';
    end
        

end

% Construct the components.
text      = uicontrol('Style','Text',...
             'String','---------------------------------------------',...
             'Position',[width/2-160 height-33 320,14],...
             'FontSize',14);
title      = uicontrol('Style','Text',...
             'String','Group Selection Tool','Position',[width/2-150 height-30 300,30],...
             'FontSize',16);


startpos = [30 height-100 bw bh];
    
%for i=1:length(Clusternames)
for i=1:num_vis
    % Clustername
    textClusternamepos = startpos;
    textClusternamepos(2) = startpos(2) - i*30;
    textClusternamepos(3) = 300;
    
    tmp = textClusternamepos;
    tmp(2) = tmp(2)-9;
    tmp(3) = 1350;
    textline      = uicontrol('Style','Text',...
        'String','_____________________________________________________________________________________________________________________________________________________________________________________________________',...
        'Position',tmp,...
        'FontSize',10);
    handles.textClustername{i}     = uicontrol('Style','edit',...
        'String',[Clusternames{i,1} ' -> (' Clustercoor{i,1} ')'],...
        'Position',textClusternamepos,...
        'FontSize',10);

end

rkCn = textClusternamepos(1)+textClusternamepos(3);
% radiobutton width
rbw = 35;


% schleife ueber die Buttons / groups
for k=1:30
    % Groupname
    groupnamepos = startpos;
    groupnamepos(1) = rkCn-10 + (k-1)*rbw;
    v = mod((k-1),2);
    groupnamepos(2) = startpos(2)-10+v*bh;
    groupnamepos(3) = 50;
    editgroupname{k}      = uicontrol('Style','edit',...
        'String','NaN',...
        'Position',groupnamepos,...
        'FontSize',6);
    handles.editgroupname{k} = editgroupname{k};
    for i=1:num_vis%length(Clusternames)
        
        
        % Radiobutton
        pos = [rkCn+10+(k-1)*rbw startpos(2) - i*30, rbw, bh];
        handles.rb{i,k}  = uicontrol('Style','radiobutton',...
            'String','|','Position',pos,...
            'FontSize',10,'Value',0);
        
        
    end
end
% 

handles.pb_saveclusterg  = uicontrol('Style','pushbutton',...
             'String','save cluster groups','Position',[20 height-20-bh bw*2 bh],...
             'FontSize',10,...
             'Callback',{@pb_saveclusterg_Callback});

handles.pb_finish  = uicontrol('Style','pushbutton',...
             'String','save cluster groups and close','Position',[40+bw*2 height-20-bh bw*2 bh],...
             'FontSize',10,...
             'TooltipString','Changed Clusternames will be saved in Cluster.mat',...
             'Callback',{@pb_finish_Callback});

handles.pb_prev  = uicontrol('Style','pushbutton',...
             'String','prev page','Position',[60+bw*4 height-20-bh bw bh],...
             'FontSize',10,...
             'TooltipString','next page',...
             'Callback',{@pb_prev_Callback});

handles.pb_next  = uicontrol('Style','pushbutton',...
             'String','next page','Position',[70+bw*5 height-20-bh bw bh],...
             'FontSize',10,...
             'TooltipString','next page',...
             'Callback',{@pb_next_Callback});



set(gcf,'CloseRequestFcn',@my_closefcn);

try
    names = Cluster{i,1}.net.global.network_names;
    for j=1:length(names)
        set(handles.editgroupname{j},'String',names{j});
    end
    
catch
end

handles.Clusternames = Clusternames;
handles.Cluster = Cluster;
handles.Clustergroupname = Clustergroupname;
handles.Clustercoor = Clustercoor;

handles.num_cluster = num_cluster;
handles.num_vis = num_vis;
handles.start_vis = 1;
handles.end_vis = num_vis;
handles.offset = 0; % offse fuer meherere Visulaisierungsseiten der Clusternamen
[ ~, handles.num_group] = size(handles.rb);
assignin('base','Clusternames',Clusternames);
%set(handles.table,'Data',[Clusternames Clustercoor Clustergroup Clustergroupname]);


% versuche die connectivitaetsmatrix zu laden
try
    handles.N_conn = handles.Cluster{1,1}.net.global.N_conn;
    %names = Cluster{1,1}.net.global.network_names;
catch
    handles.N_conn = 0;
end


guidata(f,handles)

set_radiobuttons_to_N_conn(handles);


set(f,'Name','Cluster Grouping');
set(f,'Visible','on')
% UIWAIT makes ck_rc_analyse_ICs_options wait for user response (see UIRESUME)
uiwait(f);

end

function set_radiobuttons_to_N_conn(handles)
handles = guidata(gcf);
num_group = handles.num_group;
N_conn = handles.N_conn;

% ueberprueffe die groesse
[z,s] = size(N_conn);
for i=z:handles.num_cluster
    for j=s:num_group
        N_conn(i,j)=0;
    end
end

% update die aktuell sichtbaren radiobuttons
offset = handles.offset;
num_vis = handles.num_vis;
start_vis = handles.start_vis;
end_vis = handles.end_vis;
num_cluster = handles.num_cluster;

%fprintf('update () start: %.0d   ende: %.0d    offset: %.0d \n',start_vis,end_vis,offset);
%fprintf('update () num_cluster: %.0d   num_vis: %.0d   \n',num_cluster,num_vis);

for i=start_vis:start_vis+num_vis-1
    for j=1:num_group
        if i<=end_vis
            set(handles.rb{i-offset,j},'Value',N_conn(i,j));
        else
            set(handles.rb{i-offset,j},'Value',0);
        end
    end
end

guidata(gcf,handles)
end


function pb_saveclusterg_Callback(hObject,eventdata)

handles = guidata(gcf);

% Data = get(handles.table,'Data');
if isfield(handles,'Cluster')
    Cluster = handles.Cluster;
else
    load(fullfile(handles.outdir,'Cluster.mat'));
end
num_cluster = handles.num_cluster;
num_group = handles.num_group; 

handles.networknames = get_networknames(handles);
handles.num_networknames = length(handles.networknames);
handles.N_conn = get_N_conn(handles);
assignin('base','N_conn',handles.N_conn);
% bestimme die Anzahl der Netzwerke (eintraege in der Matrix)
handles.num_global_networks = sum(sum(handles.N_conn)>0);
% vergleiche diese Zahl mit der Anzahl der eingetragenen Netzwerknamen
if handles.num_global_networks ~= handles.num_networknames
    warndlg('matrix entries did not correspond to networknames');
    return;
end

% schreibe die Infos in die Clusterstruktur
for i=1:num_cluster
    Cluster{i,1}.net.global.num_global_networks = handles.num_global_networks;
    Cluster{i,1}.net.global.network_names = handles.networknames;
    Cluster{i,1}.net.global.N_conn = handles.N_conn;
    
    index = 0;
    for idx_group = 1:num_group
%        if get(handles.rb{i,idx_group},'Value')
%        assignin('base','Clusterxx',Cluster);
        if handles.N_conn(i,idx_group)
            index = index + 1;
            Cluster{i,1}.net.local.group{index} = idx_group;
            Cluster{i,1}.net.local.groupname{index} = get(handles.editgroupname{idx_group},'String');
            % ermittelung der maximalen Anzahl an angegebenen Netzwerken
        end
    end
    Cluster{i,1}.net.local.num_local_networks = index;
end



save(fullfile(handles.outdir,'Cluster.mat'),'Cluster');
assignin('base','Clustergroup',Cluster);

handles.Cluster = Cluster;
guidata(gcf,handles);
%get(handles.editgroupname{1},'String');
%tt = get(handles.editgroupname{1},'String')
%assignin('base','tt',tt);
end


function pb_finish_Callback(hObject,eventdata)
pb_saveclusterg_Callback(hObject,eventdata);
delete(gcf);
end


function pb_grouping_Callback(hObject,eventdata)
handles = guidata(gcf);

group = str2num(get(handles.editgroup,'String'));
groupname = get(handles.editgroupname,'String');
selection = get(handles.table,'UserData');
% alte Methode ohne Batch - aber gut geschrieben 02/2012 herausgenommen

Data = get(handles.table,'Data');
%assignin('base','Data',Data);

for i=1:size(selection,1)
    Data{selection(i,1),3}=group;
    Data{selection(i,1),4}=groupname;
    
end

set(handles.table,'Data',Data);
%update(hObject,eventdata);
    
    
guidata(gcf,handles);

end


function my_closefcn(hObject,eventdata) 
delete(gcf)
end

function pb_next_Callback(hObject,eventdata)
handles = guidata(gcf);
if handles.end_vis > handles.num_cluster
    fprintf('bug in pb_next Callback end_vis to big\n');
end
if handles.end_vis == handles.num_cluster
    return;
end
if handles.end_vis < handles.num_cluster
    handles.N_conn = get_N_conn(handles);
    handles.start_vis = handles.end_vis + 1; % +1
    handles.end_vis = handles.start_vis + handles.num_vis-1;
    handles.offset= handles.start_vis - 1;
end
if handles.end_vis>handles.num_cluster
    handles.end_vis = handles.num_cluster;
end
guidata(gcf,handles);
update(handles);

end

function pb_prev_Callback(hObject,eventdata)
handles = guidata(gcf);

handles = guidata(gcf);
if handles.start_vis < 1
    fprintf('bug in pb_prevc Callback start_vis to small\n');
end
if handles.start_vis == 1
    return;
end
if handles.start_vis >1
    handles.N_conn = get_N_conn(handles);
    handles.end_vis = handles.start_vis -1;
    handles.start_vis = handles.start_vis - handles.num_vis ;
    handles.offset = handles.start_vis -1;
end

if handles.start_vis<1
     fprintf('bug in pb_prevc Callback start_vis to small\n');
end
  
guidata(gcf,handles);

update(handles);

end

% es wird die Struktur N_conn welche die Matrixeintraege fuer die 
% Radiobuttons enthaelt upgedated im Hinblicka auf die aktuellen Eintraege
% und zurueck gegeben
% die handles.N_conn Struktur wird nicht upgedated
function N_conn = get_N_conn(handles)
handles = guidata(gcf);
Cluster = handles.Cluster;
start_vis = handles.start_vis;
end_vis = handles.end_vis;
num_cluster = handles.num_cluster;
num_group = handles.num_group; 
offset = handles.offset;

if isfield(handles,'N_conn')
    N_conn = handles.N_conn;
else
    N_conn = 0;
end
for i=start_vis:end_vis
for j=1:num_group
    if get(handles.rb{i-offset,j},'Value')
        N_conn(i,j) = 1;
        %fprintf('%d  %d \n',i,j);
    else
        N_conn(i,j) = 0;
    end
end
end
assignin('base','N_conn_after_get_N_conn',N_conn);

end

% gibt eine cell Liste der Netzwerknamen zurueck
function networknames = get_networknames(handles)

handles = guidata(gcf);
num_group = handles.num_group; 
%idx = 0;
for j=1:num_group
    name = get(handles.editgroupname{j},'String');
    if strcmp(name,'NaN')
        % dann wird das Netzwerk nicht betrachtet
    else
        %idx = idx + 1;
        networknames{j} = name;
    end
end
end


function update(handles)
handles = guidata(gcf);
% update die Eintraege in den Textfeldern die die Cluster beschreiben
% pruefe ob bereits Netzwerke erstellt wurden und trage diese ein
%disp(handles.num_cluster_vis_end)
start_vis = handles.start_vis;
end_vis = handles.end_vis;
offset = handles.offset;
Cluster = handles.Cluster;
num_vis = handles.num_vis;
num_cluster = handles.num_cluster;
% schreibe die neuen Clusternamen
%size(handles.textClustername)

%fprintf('update () start: %.0d   ende: %.0d    offset: %.0d \n',start_vis,end_vis,offset);
%fprintf('update () num_cluster: %.0d   num_vis: %.0d   \n',num_cluster,num_vis);


if num_cluster<end_vis
    end_vis = num_cluster;
end


% trage die neuen Clusternamen ein
for i=start_vis:start_vis+num_vis-1 %-1 28.09.2018 ck 
    
    % Clustername
    if i<=end_vis
        handles.textClustername{i-offset}.String = handles.Clusternames{i,1};
    else
        handles.textClustername{i-offset}.String = '';
    end
end
% update die Radiobutton eintraege
N_conn = handles.N_conn;
assignin('base','N_conn',N_conn);
assignin('base','Clusterx',Cluster);
handles.networknames = get_networknames(handles);
% for j=1:size(N_conn,2)
%     %set(handles.editgroupname{j},'String',names{j});
%     for i=start_vis:end_vis%size(N_conn,1)
%         set(handles.rb{i-offset,j},'Value',N_conn(i,j));    
%     end
% end

handles.N_conn = N_conn;
handles.Cluster = Cluster;


guidata(gcf,handles);

set_radiobuttons_to_N_conn(handles);

end
