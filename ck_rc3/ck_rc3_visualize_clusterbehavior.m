function [ output_args ] = ck_rc3_visualize_clusterbehavior(outdir,G)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% die Struktur X hat die Felder

% clusternames ... die Namen der Cluster
% G1.BG behavioralen Daten der Gruppe 1
% G1.TG1 Titel der Spalten in X.BG1 als cell vector
% G2.BG behaviorale Daten der Gruppe 2
% G2.TG Titel der Spalten in X.BG2
% G1.R  R WErt der Korrelation zwischen Konnektivitaet und Behavior fuer
% eine Cluster zu Cluster verbindung
% G2.R 
% G1.P  entsprechender P Wert
% G2.P 
% G1.R_3d ... 3D Struktuer der Konnektivitaet (Area x Area x Subject)
% G2.R_3d 

% G.G1 = G1;
% G.G2 = G2;
% G.clusternames = G1.clusternames;
show(outdir,G);
 

end

function adjust_C2C_figure(clusternames)
handles = guidata(gcf);
areas=clusternames;
dim =size(areas,1);
% setze als axes die main axes
axes(handles.axmain);
set(handles.axmain,'YTick',[1:dim])
set(handles.axmain,'YTickLabel',areas)
%set(gca,'XTickLabel',{'-pi','-pi/2','0','pi/2','pi'})
set(handles.axmain,'XTick',[1:dim])
set(handles.axmain,'XTickLabel',areas')
%set(gca,'Position',[0.060,0.0600,1,0.9200]);
%set(gca,'Position',[0.05,0.05, 0.9, 0.7]);

usethesevalues=get(handles.axmain,'XTickLabel');
set(handles.axmain,'XTickLabel',[]);
putthemhere=get(handles.axmain,'XTick');
ylimits=get(handles.axmain,'Ylim');
ypos=ylimits(1)-.01*(ylimits(2)-ylimits(1));
th=text(putthemhere,ypos*ones(1,length(putthemhere)),usethesevalues,...       
'HorizontalAlignment','left','rotation',90,'parent',gca);
colorbar();
title(handles.axmain,handles.current.name1);
if handles.cb_cm.Value
    cm = eval(handles.edit_cm.String);
    %clim = caxis
    caxis(cm)
end

end


function txt = myupdatefcn2(empt,event_obj)
% Customizes text of data tips
% der listener darf nicht seine eigene axes überschreiben 
% sonst gibt es einen Fehler

handles = guidata(gcf);
%idx = handles.current.behav_idx;
pos = get(event_obj,'Position');

handles.pos = pos;
guidata(gcf,handles);
% loesche das alte REchteck
delete_old_rect();
if isfield(handles,'pos')
    % update fig 1
    show_subject_conn(handles.pos);
    % update fig 2
    txt = show_correlation(1,handles.pos);
    % update fig 3
    txt = show_correlation(2,handles.pos);
end

end




function txt = show_correlation(group,pos)
handles = guidata(gcf);
clusternames = handles.G.clusternames;     
R2_3d = handles.G.G2.R_3d;
R1_3d = handles.G.G1.R_3d;


if group==1
    axes(handles.axg1);
    subjects = [1:size(R1_3d,3)];
    handles.current.behav_idx
    b = handles.G.G1.BG(:,handles.current.behav_idx);
    c = squeeze(R1_3d(pos(2),pos(1),:));
    txt = {['Y: ',clusternames{pos(1)}],['X: ',clusternames{pos(2)}]};
    scatter(b,c,'x');
    %scatter(subjects,correll','o');
    xlabel(handles.G.G1.TG{handles.current.behav_idx});
    ylabel(handles.current.name2);
    titlestr = [ clusternames{pos(2)}, clusternames{pos(1)}, '  R = ', num2str(handles.G.G1.RR(pos(1),pos(2),handles.current.behav_idx)), ...
        '  p= ', num2str(handles.G.G1.RP(pos(1),pos(2),handles.current.behav_idx))];
    title(titlestr);
    handles.current.export.g1behav_idx = handles.current.behav_idx;
    handles.current.export.g1b = b;
    handles.current.export.g1c=c;
    handles.current.export.g1titlestr = titlestr;
    handles.current.export.g1txt = txt;
    handles.current.export.g1xlabel = handles.G.G1.TG{handles.current.behav_idx};
    handles.current.export.g1ylabel = handles.current.name2;

else
    subjects = [1:size(R2_3d,3)];
    axes(handles.axg2);
    b = handles.G.G2.BG(:,handles.current.behav_idx);
    c = squeeze(R2_3d(pos(1),pos(2),:));
    txt = {['Y: ',clusternames{pos(2)}],['X: ',clusternames{pos(1)}]};
    scatter(b,c,'x');
    %scatter(subjects,correll','o');
    handles.G.G2.TG
    xlabel(handles.G.G2.TG{handles.current.behav_idx});
    ylabel(handles.current.name2);
    titlestr = [clusternames{pos(2)}, clusternames{pos(1)}, '  R = ', num2str(handles.G.G2.RR(pos(1),pos(2),handles.current.behav_idx)), ...
        '  p= ', num2str(handles.G.G2.RP(pos(1),pos(2),handles.current.behav_idx))];
    title(titlestr);
    handles.current.export.g2behav_idx = handles.current.behav_idx;
    handles.current.export.g2b = b;
    handles.current.export.g2c = c;
    handles.current.export.g2titlestr = titlestr;
    handles.current.export.g2txt = txt;
    handles.current.export.g2xlabel = handles.G.G2.TG{handles.current.behav_idx};
    handles.current.export.g2ylabel = handles.current.name2;
end

dy = (max(c)-min(c))/40;
text(b,c-dy,cellstr(num2str(subjects')));

hold on
 
%interp1(behavioral,correll,'linear')
yp = polyfit(b,c,1);
lr = polyval(yp,b);
plot(b,lr,'r-');
lrstr = ['y = ' num2str(yp(1)) ' x + ' num2str(yp(2))];
text(min(b),polyval(yp,min(b)),lrstr);
hold off

%zeichne die Kaestchen in das imagesc
axes(handles.axmain);
if group ==1
    handles.rectg1 = rectangle('Position',[pos(2)-0.5,pos(1)-0.5,1,1],'EdgeColor','r','LineWidth',2);
    handles.liney1g1 = rectangle('Position',[pos(1)-0.5,0,1,size(handles.G.RR12,1)],'EdgeColor','r','LineWidth',0.3);
    handles.liney2g1 = rectangle('Position',[pos(2)-0.5,0,1,size(handles.G.RR12,1)],'EdgeColor','r','LineWidth',0.3);
    handles.linex1g1 = rectangle('Position',[0,pos(2)-0.5,size(handles.G.RR12,1),1.0],'EdgeColor','r','LineWidth',0.3);
    handles.linex2g1 = rectangle('Position',[0,pos(1)-0.5,size(handles.G.RR12,1),1.0],'EdgeColor','r','LineWidth',0.3);
    handles.current.export.g1yp = yp;
    handles.current.export.g1lr = lr;
    handles.current.export.g1linetext = lrstr;
else
    handles.rectg2 = rectangle('Position',[pos(2)-0.5,pos(1)-0.5,1,1],'EdgeColor','r','LineWidth',2);
    handles.liney1g2 = rectangle('Position',[pos(1)-0.5,0,1,size(handles.G.RR12,1)],'EdgeColor','r','LineWidth',0.3);
    handles.liney2g2 = rectangle('Position',[pos(2)-0.5,0,1,size(handles.G.RR12,1)],'EdgeColor','r','LineWidth',0.3);
    handles.linex1g2 = rectangle('Position',[0,pos(2)-0.5,size(handles.G.RR12,1),1.0],'EdgeColor','r','LineWidth',0.3);
    handles.linex2g2 = rectangle('Position',[0,pos(1)-0.5,size(handles.G.RR12,1),1.0],'EdgeColor','r','LineWidth',0.3);
    handles.current.export.g2yp = yp;
    handles.current.export.g2lr = lr;
    handles.current.export.g2linetext = lrstr;
end

guidata(gcf,handles);
end

function delete_old_rect()

handles = guidata(gcf);
axes(handles.axmain);
if isfield(handles,'rectg1')
    delete(handles.rectg1);
end
if isfield(handles,'linex1g1')
    delete(handles.linex1g1);
end
if isfield(handles,'liney1g1')
    delete(handles.liney1g1);
end
if isfield(handles,'linex2g1')
    delete(handles.linex2g1);
end
if isfield(handles,'liney2g1')
    delete(handles.liney2g1);
end
if isfield(handles,'rectg2')
    delete(handles.rectg2);
end
if isfield(handles,'linex1g2')
    delete(handles.linex1g2);
end
if isfield(handles,'liney1g2')
    delete(handles.liney1g2);
end
if isfield(handles,'linex2g2')
    delete(handles.linex2g2);
end
if isfield(handles,'liney2g2')
    delete(handles.liney2g2);
end

guidata(gcf,handles)

end

function pb_next_Callback(hObject,eventdata)
handles = guidata(hObject);

fprintf('next Column\n');
if handles.current.behav_idx<size(handles.G.G1.BG,2)
    handles.current.behav_idx = handles.current.behav_idx + 1;
end
guidata(hObject,handles);
update_all();
%axes(handles.axmain);
%imagesc(squeeze(handles.G.RR12(:,:,handles.current.behav_idx)));
%adjust_C2C_figure(handles.X.clusternames);
%guidata(hObject,handles);

end

function pb_export_Callback(hObject,eventdata)
handles = guidata(hObject);
exportname = handles.edit_export.String;
%T = handles.current.export;
assignin('base',handles.edit_export.String,handles.current.export);

%guidata(hObject,handles);
%update_all();
% axes(handles.axmain);
% imagesc(squeeze(handles.G.RR12(:,:,handles.current.behav_idx)));
% adjust_C2C_figure(handles.X.clusternames);
% guidata(hObject,handles);

end

function pb_show_conn_Callback(hObject,eventdata)
handles = guidata(hObject);

if ~isfield(handles,'C2C_stat')
    load(fullfile(handles.outdir,'C2C_stat.mat'));
    handles.C2C_stat = C2C_stat;

end
handles.current.showconn = 1;
handles.current.showbehav = 0;
guidata(hObject,handles);
update_all();
% axes(handles.axmain);
% handles.imagesc = imagesc(handles.C2C_stat.MeanZdiff);
% adjust_C2C_figure(handles.X.clusternames);


%guidata(hObject,handles);

end
       
function pb_show_behav_Callback(hObject,eventdata)

handles = guidata(hObject);
%show_behav();
handles.current.showconn = 0;
handles.current.showbehav = 1;
guidata(hObject,handles);
update_all();
end
% 
% function show_behav()
% handles = guidata(gcf);
% axes(handles.axmain);
% imagesc(squeeze(handles.G.RR12(:,:,handles.current.behav_idx)));
% adjust_C2C_figure(handles.X.clusternames);
% 
% 
% end

function lb_behav_Callback(hObject,eventdata)
% handles = guidata(hObject);
% handles.current.behav_idx = get(hObject,'Value');
% guidata(hObject,handles);
% %show_behav();
% if isfield(handles,'pos')
% txt = show_correlation(2,handles.pos);
% txt = show_correlation(1,handles.pos);
% show_subject_conn(handles.pos);
% end
update_all();
end


function lb_conn_Callback(hObject,eventdata)
%handles = guidata(hObject);
% handles.conn_idx = get(hObject,'Value');
% guidata(hObject,handles);

%adjust_C2C_figure(handles.X.clusternames);
%guidata(hObject,handles);
update_all();

end


function pb_show_power_Callback(hObject,eventdata)
handles = guidata(hObject);
%estimate power
for i=1:length(handles.Cluster)
    %p(i)=Cluster{i}.
end
handles.axpower = axes('Position',[0.70 0.05 0.05 0.85] );
axes(handel.axpower);

if isfield(handles,'pos')
txt = show_correlation(2,handles.pos);
txt = show_correlation(1,handles.pos);
show_subject_conn(handles.pos);
end
end


function show(outdir,G)

fprintf('start');
handles.outdir = outdir;
handles.G = G;
screensize = get( 0, 'Screensize' );
screensize(4) = screensize(4)-30;
handles.current.behav_idx = 1; % der index der behavioralen Daten der gerade gzeigt wird

f = figure('Visible','off','Position',screensize,...
    'MenuBar','none',...
    'Resize','off',...
    'Toolbar','none');

handles.axg2 = axes('Position',[0.78 0.15 0.2 0.24 ]);
handles.axg1 = axes('Position',[0.78 0.45 0.2 0.24 ]);
handles.axconn = axes('Position',[0.78 0.76 0.2 0.22 ]);


%set(handles.imagesc,'ButtonDownFcn',@mydatacursorupdate);

%set(dcm_obj.CurrentDataCursor, 'Marker','o', 'MarkerFaceColor','b');

% handles.pb_next = uicontrol('Style','pushbutton',...
%               'String','next','Position',[screensize(1)-100 screensize(2)-50 40 30],...
%               'FontSize',10,...
%               'Callback',{@pb_next_Callback});
bw = 40;
bh = 25;
lbw = 130;
lbh = 130;

% listbox mit behavioralen Daten
handles.lb_behav = uicontrol('Style','listbox',...
              'String',handles.G.G1.TG,...
              'Position',[f.Position(3)-200 20 lbw lbh],...
              'FontSize',10,...
              'Callback',{@lb_behav_Callback});

handles.show_var_list = {'conn vs. behav'; 'brain conn'};
handles.lb_show = uicontrol('Style','listbox',...
              'String',handles.show_var_list,...
              'Position',...
              [handles.lb_behav.Position(1)-lbw-10,...
              handles.lb_behav.Position(2),lbw,lbh],...
              'FontSize',10,...
              'Callback',{@lb_conn_Callback});

% Listbox mit conn Daten
handles.connvarlist = get_C2C_stat_var_list_to_display(handles.outdir);

handles.lb_conn = uicontrol('Style','listbox',...
              'String',handles.connvarlist,...
              'Position',...
              [handles.lb_show.Position(1)-lbw-10,...
              handles.lb_show.Position(2),lbw,lbh],...
              'FontSize',10,...
              'Callback',{@lb_conn_Callback});
          
handles.edit_export = uicontrol('Style','edit',...
              'String','[-0.5 0.5]',...
              'Position',...
              [handles.lb_behav.Position(1),...
              handles.lb_behav.Position(2)+10+lbh,bw*2,bh],...
              'FontSize',10,...
              'TooltipString','Variable name for Export to workspace'); 
 handles.pb_export = uicontrol('Style','pushbutton',...
               'String','export','Position',[handles.edit_export.Position(1)+10+bw*2,...
               handles.edit_export.Position(2),bw, bh],...
               'FontSize',10,...
               'TooltipString','export to workspace',...
               'Callback',{@pb_export_Callback});

          

handles.edit_cm = uicontrol('Style','edit',...
              'String','[-0.5 0.5]',...
              'Position',...
              [handles.lb_conn.Position(1)-lbw-10,...
              handles.lb_conn.Position(2),bw*2,bh],...
              'FontSize',10,...
              'TooltipString','restrict the colormap between values');   

handles.cb_cm = uicontrol('Style','checkbox',...
              'String','restrict colormap',...
              'Position',...
              [handles.edit_cm.Position(1),...
              handles.edit_cm.Position(2)+bh+10,bw*3,bh],...
              'FontSize',10,...
              'TooltipString','restrict the colormap between values');   
 handles.pb_visualizeCluster = uicontrol('Style','pushbutton',...
               'String','vis','Position',[handles.cb_cm.Position(1),...
               handles.cb_cm.Position(2)-2*bh,bw, bh],...
               'FontSize',10,...
               'TooltipString','visualize Cluster',...
               'Callback',{@pb_visualizeCluster_Callback});
          
  

         
          
table_pos = [handles.axg1.Position(1),...
    handles.axg1.Position(2)+handles.axg1.Position(4)+50,...
    handles.axg1.Position(3),...
    f.Position(4)-handles.axg1.Position(2)-handles.axg1.Position(4)-100];
 
% table = uitable('Parent', f, 'Data', cell(length(N.networknumbers),8),...
%      'Position', table_pos ,...
%      'ColumnName',{'network name','num','p-value','p G1<G2','p G1>G2','MeanZdiff', 'Zstd', 'MeanRdiff'},...
%      'ColumnEditable',true,...
%      'ColumnWidth',{150 50 70 70 70 70 70 70},...
%      'CellSelectionCallback',@(src,evnt)set(src,'UserData',evnt.Indices));




%handles.G = G;
% define Axis
%assignin('base','RR12',G.RR12);
handles.axmain = axes('Position',[0.05 0.05 0.7 0.85]);
handles.imagesc = imagesc(squeeze(handles.G.RR12(:,:,handles.current.behav_idx)));
handles.current.name1 = 'Pearson Mean R Diff';
handles.current.name2 = 'Pearson r-value';
handles.current.imageMatrix = handles.G.RR12;

set(handles.axmain,'Units','pixels');

guidata(f,handles)
adjust_C2C_figure(G.clusternames);
dcm_obj = datacursormode(f);
set(dcm_obj,'DisplayStyle','datatip','Enable','on','UpdateFcn',@myupdatefcn2)


set(f,'Name','Behavioral correlations');

set(f,'Visible','on')
% UIWAIT makes ck_rc_analyse_ICs_options wait for user response (see UIRESUME)
%uiwait(f);


end



function my_closefcn(hObject,eventdata) 
delete(gcf)
end



function show_subject_conn(pos)

handles = guidata(gcf);
G = handles.G;
name = [G.clusternames{pos(1)} ' vs. ' G.clusternames{pos(2)}];

%'show subject conn'

c1 = 'red'; %color 1
c2 = 'blue'; %color 2
axes(handles.axconn);
cla(handles.axconn);
%figure('Name',netname);

R1 = G.G1.R_3d;
R2 = G.G2.R_3d;

if ~isfield(handles,'C2C_stat_g1')
    load(fullfile(handles.outdir,'C2C_stat_g1.mat'));
    handles.C2C_stat_g1 = C2C_stat_g1;
end
if ~isfield(handles,'C2C_stat_g2')
    load(fullfile(handles.outdir,'C2C_stat_g2.mat'));
    handles.C2C_stat_g2 = C2C_stat_g2;
end
if ~isfield(handles,'C2C_stat')
    load(fullfile(handles.outdir,'C2C_stat.mat'));
    handles.C2C_stat = C2C_stat;
end


% Personen Gruppe 1
y1=squeeze(R1(pos(1),pos(2),:));
%N.RGN1(:,num);
x1=[1:length(y1)];
%e1= N.RGN1std(:,num);
e1= y1*0;
e1n= e1*00.1; % setze die Errorbar nach unten auf Null


%Personen Gruppe 2
y2=squeeze(R2(pos(1),pos(2),:));
x2=[1:length(y2)]+length(y1)+3;
e2= y2*0;
e2n= e1*00.1; % setze die Errorbar nach unten auf Null

[p1 t] = barwitherr(e1,x1,y1,c1); %,'DisplayName',{'group1','e1'});


hold on
p2=barwitherr(e2,x2,y2,c2); %,'DisplayName','group1');

%legend([p1 t p2],'group1','std within Network','group2');

hold on 
xl1 = [x1(1):1:x2(end)];
yl1 = ones(length(xl1))*mean(squeeze(R1(pos(1),pos(2),:)));
% 
% assignin('base','y1',y1);
% assignin('base','x1',x1);
% assignin('base','e1',e1);
% assignin('base','y2',y2);
% assignin('base','x2',x2);
% assignin('base','e2',e2);
% assignin('base','xl1',xl1);
% assignin('base','yl1',yl1);
plot(xl1,yl1,c1);

xl2 = [x1(1):1:x2(end)];
yl2 = ones(length(xl1))*mean(squeeze(R2(pos(1),pos(2),:)));
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
    elseif i<=length(x1)+3
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
str = sprintf('R1 mean = %.4f  R2 mean = %.4f  R1-R2= %.4f    p = %.7f',...
    handles.C2C_stat.R1mean(pos(1),pos(2)),...
    handles.C2C_stat.R2mean(pos(1),pos(2)),...
    handles.C2C_stat.R1mean(pos(1),pos(2))-handles.C2C_stat.R2mean(pos(1),pos(2)),...
    handles.C2C_stat.Zp(pos(1),pos(2)));
xlabel(str);
ylabel('R-value');
str = sprintf('connectivity between areas: %s',name);
title(str,'Color','red','FontSize',12); 
assignin('base','y1',y1);
assignin('base','y2',y2);

end

function varlist = get_C2C_stat_var_list_to_display(outdir)

load(fullfile(outdir,'C2C_stat.mat'));
varlist = cell(0);
if isfield(C2C_stat,'MeanZdiff')
    varlist{end+1} = 'MeanZdiff';
end
if isfield(C2C_stat,'MeanSyndiff')
    varlist{end+1} = 'MeanSyndiff';
end
if isfield(C2C_stat,'MeanPLVdiff')
    varlist{end+1} = 'MeanPLVdiff';
end
end



function update_all(varargin)

handles = guidata(gcf);

% falls der string 'not main' uebergeben wird dann nicht das main fig
% ueberschreiben
if nargin==1
    str = varargin{1};
else
    str = 'all';
end
% get infos
handles.current.behav_idx = handles.lb_behav.Value;
handles.current.behav_idx
handles.current.show_idx = handles.lb_show.Value; % 1 ist behavior, 2 ist brain con
handles.current.conn_idx = handles.lb_conn.Value; % 1 ist behavior, 2 ist brain con
connlist = handles.lb_conn.String;
handles.current.conn_string = connlist{handles.current.conn_idx};


% update main image matrix 
%guidata(gcf,handles);
handles.current = get_main_image_matrix(handles);

%handles = guidata(gcf);


assignin('base','handles',handles);
if strcmp(str,'not main')
    fprintf('main window will not be updated');
else
    %axes(handles.axmain);
    imagesc(handles.current.M);
    adjust_C2C_figure(handles.G.clusternames);
end
guidata(gcf,handles);
'after imagesc'

if isfield(handles,'pos')
    'update corr'
% update fig 1
show_subject_conn(handles.pos);

% update fig 2 
txt = show_correlation(1,handles.pos);

% update fig 3
txt = show_correlation(2,handles.pos);
end

end


function current = get_main_image_matrix(handles)

%handles = guidata(gcf);
if ~isfield(handles,'C2C_stat')
    load(fullfile(handles.outdir,'C2C_stat.mat'));
    handles.C2C_stat = C2C_stat;
end

% nur brain connectivity
if handles.current.show_idx ==2 
    switch handles.current.conn_string
    case 'MeanZdiff'
%        handles.imagesc = imagesc(handles.C2C_stat.MeanZdiff);
        handles.current.name1 = 'Pearson Mean Z Diff';
        handles.current.name2 = 'Pearson r-value';
        handles.current.M = handles.C2C_stat.MeanZdiff;
        %handles.current.imageMatrixConn = handles.G.Mean;
    case 'MeanSyndiff'
%        handles.imagesc = imagesc(handles.C2C_stat.MeanSyndiff);
        handles.current.name1 = 'Synchronization';
        handles.current.name2 = 'Kopplungsindex';
        handles.current.M = handles.C2C_stat.MeanSyndiff;
    case 'MeanPLVdiff'
%        handles.imagesc = imagesc(handles.C2C_stat.MeanPLVdiff);
        handles.current.name1 = 'Phase Locking Value';
        handles.current.name2 = 'Phase Locking Value';
        handles.current.M = handles.C2C_stat.MeanPLVdiff;
    otherwise
    end
end

if handles.current.show_idx ==1 
% hier das ausgewaehlte image rein
switch handles.current.conn_string
    case 'MeanZdiff'
%        handles.imagesc = imagesc(handles.C2C_stat.MeanZdiff);
        handles.current.name1 = 'Pearson Mean Z Diff';
        handles.current.name2 = 'Pearson r-value';
        handles.current.M = squeeze(handles.G.RR12(:,:,handles.current.behav_idx));
        %handles.current.imageMatrixConn = handles.G.Mean;
    case 'MeanSyndiff'
%        handles.imagesc = imagesc(handles.C2C_stat.MeanSyndiff);
        handles.current.name1 = 'Synchronization';
        handles.current.name2 = 'Kopplungsindex';
        handles.current.M = squeeze(handles.G.SynR12(:,:,handles.current.behav_idx));
    case 'MeanPLVdiff'
%        handles.imagesc = imagesc(handles.C2C_stat.MeanPLVdiff);
        handles.current.name1 = 'Phase Locking Value';
        handles.current.name2 = 'Phase Locking Value';
        handles.current.M = squeeze(handles.G.PLVR12(:,:,handles.current.behav_idx));
    otherwise
end


end
%guidata(gcf,handles);
current = handles.current;
end





function pb_visualizeCluster_Callback(hObject, eventdata)
'visualize' 
handles = guidata(hObject);
if ~isfield(handles,'Cluster')
    load(fullfile(handles.outdir,'Cluster.mat'));
    handles.Cluster = Cluster;
end

% if size(handles.Cluster,1) == 0
%     fprintf('No clusters loaded!\n');
%     return
% end
    handles.pos(1)
    handles.pos(2)
    %Cluster{handles.pos(1)}.name
aktdir = cd;
    if ~isfield(handles,'SPM_visualize')
       load(which('sb_roi_cluster_SPM_template.mat'));

       if  exist('SPM')
           handles.SPM_visualize.SPM = SPM;
           handles.SPM_visualize.VOL = VOL;
           fprintf('loaded: %s\n',which('sb_roi_cluster_SPM_template.mat'));
           
       else
           fprintf('something went wrong reading: %s\n',which('sb_roi_cluster_SPM_template.mat'));
       end

    end
    assignin('base','handles',handles);
    handles.pos(1)
    handles.pos(2)
    %val = get(handles.listbox_clusterlist,'Value');
    %str = get(handles.listbox_clusterlist,'String');
    %val = val(1); %no multiselection supported here!
    
    val = handles.pos(1);
    str = handles.Cluster{val}.name;
    
    SPM = handles.SPM_visualize.SPM;
    %SPM.swd = handles.working_dir ;
    SPM.swd = handles.outdir ;
    VOL = handles.SPM_visualize.VOL;
    %VOL.swd = handles.working_dir ;

    Cluster = handles.Cluster{val,1};
%    Cluster = handles.Cluster{handles.pos(1),1};
 if size(Cluster.Z,2)   == 0
     errordlg ('Cluster is empty!');
     return
 end
    VOL.Z = Cluster.Z;
    VOL.XYZmm = Cluster.XYZ;
    VOL.XYZ = Cluster.XYZ;

    VOL.mat = VOL.M;
    VOL.title = str;
    VOL.STATstr = sprintf('%i clusters\nmin-z:%g\nmax-z:%g\nmean-z:%g',size(Cluster.XYZ,2),min(VOL.Z),max(VOL.Z),mean(VOL.Z));
    VOL.u = min(VOL.Z);
    
    for i = 1:size(Cluster.XYZ,2)
        j = Cluster.XYZ(1,i);
        k = Cluster.XYZ(2,i);
        l = Cluster.XYZ(3,i);
        [x,y,z] = mni2nifti(VOL.mat,j,k,l);
        if x>0 && y > 0 && z > 0
           % Y(x,y,z) = Cluster.Z(i);
            VOL.XYZ(1,i) = x;
            VOL.XYZ(2,i) = y;
            VOL.XYZ(3,i) = z;
        end

    end
    
    %add additional point at 1x1x1 with doubled z-value for legend-spreading ;-)
    if min(VOL.Z) == max(VOL.Z)
        VOL.XYZ(1,i+1) = 1;
        VOL.XYZ(2,i+1) = 1;
        VOL.XYZ(3,i+1) = 1;
        VOL.Z(1,i+1)  = 2*VOL.Z(1,1);
        [x,y,z] = nifti2mni(VOL.mat,1,1,1);
            VOL.XYZmm(1,i+1) = x;
            VOL.XYZmm(2,i+1) = y;
            VOL.XYZmm(3,i+1) = z;
    end
        
    %%%%%%%%%%%%%%%%%
    
    [hReg,xSPM,Finter,Fgraph] = sb_spm_results_ui_new(SPM,VOL);
  try
    spm_sections(xSPM,hReg,handles.t1_file);
    assignin('base','Cluster',Cluster)
    pos = find(Cluster.Z == max(Cluster.Z));
    last_MNI = [Cluster.XYZ(1,pos(1)) Cluster.XYZ(2,pos(1)) Cluster.XYZ(3,pos(1))];
    spm_orthviews('SetCoords',last_MNI,hReg);
	spm_orthviews('Reposition');
    handles.last_MNI = last_MNI;
  catch
     fprintf('something went wrong 250620150945\n'); 
  end
    handles.cf = gcf;
    handles.output = hObject;
    guidata(hObject, handles);
    cd(aktdir);
end


function [x y z] = mni2nifti(mat, x, y, z)

xyz = round( [x y z 1] * inv(mat') );
x = xyz(1);
y = xyz(2);
z = xyz(3);
end

function [x y z] = nifti2mni(mat, x, y, z)

xyz = round( [x y z 1] * (mat') );
x = xyz(1);
y = xyz(2);
z = xyz(3);
end