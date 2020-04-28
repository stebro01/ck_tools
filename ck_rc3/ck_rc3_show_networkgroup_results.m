function [ output_args ] = ck_rc3_show_networkgroup_results(outdir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% figure mit table
% name des Netzwerks, nummer, Mittelwert gruppe1, gruppe2, differenz,
% std1,std2, p Wert, signifikazen

show_networkgroup_table(outdir);
 

end



function show_networkgroup_table(outdir)

width = 800;
height = 1000;
f = figure('Visible','off','Position',[50,50,800,height],...
    'MenuBar','none',...
    'Resize','off',...
    'Toolbar','none');



%Define Buttons
bw = 90;
bh = 20;
%Define Editfields
ew = 80;
eh = 20;
% GUI datan
handles = struct;
handles.output = struct;
handles.output.cancel = 0;
handles.outdir = outdir;


try
    load(fullfile(handles.outdir,'C2C_stat.mat'));
    load(fullfile(handles.outdir,'IntraNetwork.mat'));
catch 
    warndlg('no C2C_stat.mat file available in current directory','Warning');
    handles.output.cancel = 1;
    guidata(hObject,handles);
    return;
end
handles.C2C_stat = C2C_stat.C2C;
% if ~isfield(C2C_stat,'Network')
%     warndlg('you have to specify networks first');
%     return;
% end
N = IntraNetwork;
handles.N = N;
%assignin('base','Network',N);

% 
% % converting Cluster
% for i=1:length(Cluster)
%     Clusternames{i,1} = Cluster{i}.name;
%     middle = floor(length(Cluster{i}.XYZ)/2);
%     Clustercoor{i,1} = num2str(round(Cluster{i}.XYZ(:,middle)'));
%     if isfield(Cluster{i},'group')
%         Clustergroup{i,1} = Cluster{i}.group;
%     else
%         Clustergroup{i,1} = 1;
%     end
% end


str = sprintf('Network Results');
% Construct the components.
text      = uicontrol('Style','Text',...
             'String','---------------------------------------------',...
             'Position',[width/2-160 height-33 320,14],...
             'FontSize',14);
title      = uicontrol('Style','Text',...
             'String',str,'Position',[width/2-150 height-30 300,30],...
             'FontSize',16);


startpos =  [40,220,ew,eh];
dy = 50;

% textepipos = startpos;
% textepi      = uicontrol('Style','Text',...
%              'String','EPI directory',...
%              'Position',textepipos,...
%              'FontSize',8);


% f = figure('Position', [100 100 752 350]);

table_pos = [50 height-70-(length(N.Zp)+1)*23 680 (length(N.Zp)+1)*23 ];
table2_pos = [50 height-70-280 680 280];

table = uitable('Parent', f, 'Data', cell(length(N.networknumbers),8),...
    'Position', table_pos ,...
    'ColumnName',{'network name','num','p-value','p G1<G2','p G1>G2','MeanZdiff', 'Zstd', 'MeanRdiff'},...
    'ColumnEditable',true,...
    'ColumnWidth',{150 50 70 70 70 70 70 70},...
    'CellSelectionCallback',@(src,evnt)set(src,'UserData',evnt.Indices));

table2_pos = [50 height-70-(length(N.Zp)+1)*23 680 (length(N.Zp)+1)*23 ];
table2_pos = [50 height-70-350-560 680 560];

table2 = uitable('Parent', f, 'Data', cell(length(N.networknumbers),6),...
    'Position', table2_pos,...
    'ColumnName',{'network name','group','R-mean','R-std','Z-mean','Z-std',},...
    'ColumnEditable',true,...
    'ColumnWidth',{150 50 70 70 70 70 },...
    'CellSelectionCallback',@(src,evnt)set(src,'UserData',evnt.Indices));



% 
pb_show_subjects  = uicontrol('Style','pushbutton',...
              'String','show subjects for one net','Position',[10 height-30 bw*2 bh],...
              'FontSize',10,...
              'Callback',{@pb_show_subjects_Callback});

pb_show_subjects_all  = uicontrol('Style','pushbutton',...
              'String','show subjects for all net','Position',[10 height-60 bw*2 bh],...
              'FontSize',10,...
              'Callback',{@pb_show_subjects_all_Callback});

pb_show_network_connections  = uicontrol('Style','pushbutton',...
              'String','show network connections','Position',[width-bw*2-10 height-30 bw*2 bh],...
              'FontSize',10,...
              'Callback',{@pb_show_network_connections_Callback});

          % 
% pb_saveclusternames  = uicontrol('Style','pushbutton',...
%              'String','change cluster names','Position',[50+40+3*bw 50 bw*2 bh],...
%              'FontSize',10,...
%              'TooltipString','Changed Clusternames will be saved in Cluster.mat',...
%              'Callback',{@pb_saveclusternames_Callback});
% 

set(gcf,'CloseRequestFcn',@my_closefcn);






% Update handles structure
handles.table = table;
handles.table2 = table2;
%handles.editgroup = editgroup;

%assignin('base','Clusternames',Clusternames);
set(handles.table,'Data',[N.networknames' num2cell(N.networknumbers)',...
    num2cell(N.Zp)', num2cell(N.Z1k2p)', num2cell(N.Z2k1p)',...
    num2cell(N.MeanZdiff)', num2cell(N.Zsd)', num2cell(N.MeanRdiff)']);

g = repmat([1 2],1,length(N.Zp));
for i=1:length(N.Zp)
    names{1,i*2-1}=N.networknames{1,i};
    names{1,i*2}=N.networknames{1,i};
    Rmean{1,i*2-1}=mean(N.RGN1(:,i));
    Rmean{1,i*2}=mean(N.RGN2(:,i));
    Rstd{1,i*2-1}=std(N.RGN1(:,i));
    Rstd{1,i*2}=std(N.RGN2(:,i));
    Zmean{1,i*2-1}=mean(N.ZGN1(:,i));
    Zmean{1,i*2}=mean(N.ZGN2(:,i));
    Zstd{1,i*2-1}=std(N.ZGN1(:,i));
    Zstd{1,i*2}=std(N.ZGN2(:,i));
    
end

set(handles.table2,'Data',[names' num2cell(g)',...
    Rmean', Rstd', Zmean',Zstd']);

guidata(f,handles)


set(f,'Name','Network results');
set(f,'Visible','on')
% UIWAIT makes ck_rc_analyse_ICs_options wait for user response (see UIRESUME)
%uiwait(f);


end


function pb_show_subjects_all_Callback(hObject,eventdata)
handles = guidata(gcf);

show_network_subjects_all(handles.N);

guidata(gcf,handles);
end



function pb_show_subjects_Callback(hObject,eventdata)
handles = guidata(gcf);

%group = str2num(get(handles.editgroup,'String'));
%groupname = get(handles.editgroupname,'String');
selection = get(handles.table,'UserData');
% alte Methode ohne Batch - aber gut geschrieben 02/2012 herausgenommen

Data = get(handles.table,'Data');
%assignin('base','Data',Data);

if size(selection,1)==0
    warndlg('you have to select an network');
    return;
end
netname = Data{selection(1,1),1};
num  = Data{selection(1,1),2};

%[D] = get_network_subject_data(handles.N,name,num);
show_network_subjects(handles.N,num,netname);

guidata(gcf,handles);
end

function pb_show_network_connections_Callback(hObject,eventdata)
handles = guidata(hObject);
%group = str2num(get(handles.editgroup,'String'));
%groupname = get(handles.editgroupname,'String');
selection = get(handles.table,'UserData');
% alte Methode ohne Batch - aber gut geschrieben 02/2012 herausgenommen

Tabdata = get(handles.table,'Data');
%assignin('base','Data',Data);

if size(selection,1)==0
    warndlg('you have to select an network');
    return;
end
netname = Tabdata{selection(1,1),1};
num  = Tabdata{selection(1,1),2};

%[D] = get_network_subject_data(handles.N,name,num);
show_network_connections(handles.N,num,netname,handles.C2C_stat);



guidata(hObject,handles);
end


function show_network_connections(N,num,netname,C2C_stat)


assignin('base','N',N);
assignin('base','num',num);
assignin('base','netname',netname);
assignin('base','C2C_statn',C2C_stat);
c1 = 'red'; %color 1
c2 = 'blue'; %color 2

ni = N.network_indices{num};
nnames = N.networknames;
R1mean = C2C_stat.R1mean; 
R1std = C2C_stat.R1std; 
R2mean = C2C_stat.R2mean; 
R2std = C2C_stat.R2std; 
Zp = C2C_stat.Zp; 

r1mni = R1mean(ni,ni)
r1sni = R1std(ni,ni)
r2mni = R2mean(ni,ni)
r2sni = R2std(ni,ni)
%r2sni = ones(2,2)*0.3
zp12  = Zp(ni,ni);
%%%%%%%%%%%%
w = 800;
h = 1050;
f=figure('Name',netname,'Position',[10 10 w h]);
handles = guidata(f);
handles.C2C_stat = C2C_stat;
handles.ni = ni;
% Construct the components.
str = sprintf('Results of the Network: %s',netname);
text      = uicontrol('Style','Text',...
             'String','---------------------------------------------',...
             'Position',[w/2-160 h-33 320,14],...
             'FontSize',14);
title      = uicontrol('Style','Text',...
             'String',str,'Position',[w/2-150 h-30 300,30],...
             'FontSize',16);

% konstruiere Daten fuer die Tabelle
Tab = cell(size(r1mni,1));
Colnames = N.networknames;

%for i=1:
    
tablex_pos = [10 h-250-80 w-20 250];
tablex = uitable('Parent', f, 'Data', cell(size(r1mni)),...
    'Position', tablex_pos,...
    'ColumnName',C2C_stat.clusternames(ni)',...
    'RowName',C2C_stat.clusternames(ni)',...
    'CellSelectionCallback',{@uitable_CellSelection_Callback});

   % 'CellSelectionCallback',@(src,evnt)set(src,'UserData',evnt.Indices));

handles.tablex = tablex;
%set(handles.tablex,'Data',[num2cell(r1mni)]);
set(handles.tablex,'Data',[num2cell(zp12)]);
assignin('base','Network_p',zp12);
%h1_pos = [tablex_pos(1), tablex_pos(2)-tablex_pos(4)-30, w-20, 250];
h1_pos = [10 100 300 200];
%ax1 = axes('Position',h1_pos,'Units','pixel');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Darstellung jeder Verbindung des Netzwerks fuer beide Gruppen als
% vergleich in der axis ax1
idx = 1;
gruppe = 1;
for i=1:size(r1mni,1)
    for j=i+1:size(r1mni,2)
        y(idx,1) = r1mni(i,j); 
        y(idx,2) = r2mni(i,j); 
        e(idx,1) = r1sni(i,j); 
        e(idx,2) = r2sni(i,j); 
        %nameconn{idx} = [C2C_stat.clusternames{ni(i)} ' vs. ' C2C_stat.clusternames{ni(j)}];
        nameconn{idx} = [num2str(i) ' vs. ' num2str(j)];
        idx = idx +1;
    end
end
assignin('base','y',y);
assignin('base','e',e);
%f1 = figure;

%%%%
% problem mit den Farben wenn nur 2 Areale ein netzwerk bilden
% workaround
if size(y,1)==1 && size(y,2)==2
    y(2,1)=0;
    y(2,2)=0;
    e(2,1)=0;
    e(2,2)=0;
end

subplot(3,1,2);
[b] = barwitherr(e,y);
set(b(1),'FaceColor','red');
set(b(2),'FaceColor','blue');
set(gca,'xticklabels',nameconn);
%legend([p1 t p2],'group1','std within Network','group2');

%text('left group 1, right group2, error bars mark the within network standard deviation');
xlabel('connections ...  error bars mark the group standard deviation');
ylabel('R-value');
%str = sprintf('within network connectivity of network: ');
%title(str,'Color','red','FontSize',12); 
guidata(f,handles);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function uitable_CellSelection_Callback(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Darstellung einer Verbindung getrennt fuer alle Probanden
% als 
% h2_pos = [h1_pos(1), h1_pos(2)-h1_pos(4)-30, w-20, 250];
% %ax2 = axes('Position',h2_pos);
c1 = 'red'; %color 1
c2 = 'blue'; %color 2
handles = guidata(gcf);
fprintf('Selection\n')
pos = eventdata.Indices;
C2C_stat = handles.C2C_stat;
% %get(handles.tablex,'UserData',eventdata.Indices)
% %get(handles.tablex,'Data')
% % t = handles.tablex;
% % tablex.eventdata
% get(handles.tablex,'Data');
% selection = get(handles.tablex,'UserData')
% alte Methode ohne Batch - aber gut geschrieben 02/2012 herausgenommen
if size(pos,1)==0
    warndlg('you have to select an network');
    return;
end
posx = handles.ni(pos(1));
posy = handles.ni(pos(2));

%Tabdata = get(handles.tablex,'Data');
%assignin('base','selection',selection);


%netname = Tabdata{selection(1,1),1};
%num  = Tabdata{selection(1,1),2};
%[posx, posy] = selection(1,1);

ri1 = squeeze(C2C_stat.data.R1(posx,posy,:));
ri2 = squeeze(C2C_stat.data.R2(posx,posy,:));

assignin('base','ri1',ri2);

%figure;
subplot(3,1,3);
bar([ri1' 0 0 0 ri2']);
% bar(y);

xlabel('individuals');
ylabel('R-value');


% Personen Gruppe 1
y1=ri1';
x1=[1:length(y1)];

%Personen Gruppe 2
y2=ri2';
x2=[1:length(y2)]+length(y1)+3;

bar(x1,y1,c1); %,'DisplayName',{'group1','e1'});

hold on
bar(x2,y2,c2); %,'DisplayName','group1');

%legend([p1 t p2],'group1','std within Network','group2');

hold on 
xl1 = [x1(1):1:x2(end)];
yl1 = ones(length(xl1))*mean(ri1);
plot(xl1,yl1,c1);

xl2 = [x1(1):1:x2(end)];
yl2 = ones(length(xl1))*mean(ri2);
plot(xl2,yl2,c2);

%text('left group 1, right group2, error bars mark the within network standard deviation');
xlabel('subjects for a specific connection');
ylabel('R-value');
str = sprintf('single connection for all subjects: ');
title(str,'Color','red','FontSize',12); 
hold off



end

function show_network_subjects(N,num,netname)

c1 = 'red'; %color 1
c2 = 'blue'; %color 2

figure('Name',netname);

% Personen Gruppe 1
y1=N.RGN1(:,num);
x1=[1:length(y1)];
e1= N.RGN1std(:,num);
e1n= e1*00.1; % setze die Errorbar nach unten auf Null


%Personen Gruppe 2
y2=N.RGN2(:,num);
x2=[1:length(y2)]+length(y1)+3;
e2= N.RGN2std(:,num);
e2n= e2*0.01;

[p1 t] = barwitherr(e1,x1,y1,c1); %,'DisplayName',{'group1','e1'});


hold on
p2=barwitherr(e2,x2,y2,c2); %,'DisplayName','group1');

%legend([p1 t p2],'group1','std within Network','group2');

hold on 
xl1 = [x1(1):1:x2(end)];
yl1 = ones(length(xl1))*mean(N.RGN1(:,num));
plot(xl1,yl1,c1);

xl2 = [x1(1):1:x2(end)];
yl2 = ones(length(xl1))*mean(N.RGN2(:,num));
plot(xl2,yl2,c2);


%%% Beschriftung der einzelnen Balken durch Zahlen
for i=1:(length(x1)+length(x2)+3)
    if i<=length(x1)
        j=i;
        if y1(j)>=0
            yh = -0.05;
        else
            yh = 0.05;
            
        end
    elseif i<=length(x1)+3;
        j=0;
        yh = 0;
    else
        j=i-length(x1)-3;
        if y2(j)>=0
            yh = -0.05;
        else
            yh = 0.05;
        end
    end
    
        
    
    strx = sprintf('%s',num2str(j));
    if j==0
        strx=' ';
    end
    if mod(i,2)==0
        yh = yh *2;
    else
%        yh = 0.05;
    end
    if j<10
        xi=i-0.25;
    else
        xi = i-0.45;
    end
    
    text(xi,yh,strx,'Color','black','FontSize',14,'FontWeight','bold');
end

%text('left group 1, right group2, error bars mark the within network standard deviation');
xlabel('subjects   error bars mark the within network standard deviation');
ylabel('R-value');
str = sprintf('within network connectivity of network: %s',netname);
title(str,'Color','red','FontSize',12); 
assignin('base','y1',y1);
assignin('base','y2',y2);

end



function show_network_subjects_all(N)

netnames = N.networknames;
net_num = length(netnames);
figure('Name','Netzwerke');
for num=1:net_num
    netname = netnames{num};
subplot(net_num,1,num);

    c1 = 'red'; %color 1
c2 = 'blue'; %color 2


% Personen Gruppe 1
y1=N.RGN1(:,num);
x1=[1:length(y1)];
e1= N.RGN1std(:,num);
e1n= e1*00.1; % setze die Errorbar nach unten auf Null


%Personen Gruppe 2
y2=N.RGN2(:,num);
x2=[1:length(y2)]+length(y1)+3;
e2= N.RGN2std(:,num);
e2n= e2*0.01;

[p1 t] = barwitherr(e1,x1,y1,c1); %,'DisplayName',{'group1','e1'});


hold on
p2=barwitherr(e2,x2,y2,c2); %,'DisplayName','group1');

%legend([p1 t p2],'group1','std within Network','group2');

hold on
xl1 = [x1(1):1:x2(end)];
yl1 = ones(length(xl1))*mean(N.RGN1(:,num));
plot(xl1,yl1,c1);

xl2 = [x1(1):1:x2(end)];
yl2 = ones(length(xl1))*mean(N.RGN2(:,num));
plot(xl2,yl2,c2);


for ii=1:5:x2(end)
    if ii<x2(1)
        xnames{ii} = num2str(ii);
    else
        xnames{ii} = num2str(ii-x2(1)+1);
    end
end
ax = gca;
ax.XTick = 1:1:x2(end);
ax.XTickLabel = xnames;
%set(gca,'xticklabels',xnames);


%text('left group 1, right group2, error bars mark the within network standard deviation');
xlabel('subjects   error bars mark the within network standard deviation');
ylabel('R-value');
str = sprintf('within network connectivity of network: %s',netname);
title(str,'Color','red','FontSize',12); 
hold off
end
end


% get the subject specific data from one Network
function D = get_network_subject_data(N,name,num)

% indiviual Rmean, Rstd    von beiden Gruppen


end


function tmp()
        
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

