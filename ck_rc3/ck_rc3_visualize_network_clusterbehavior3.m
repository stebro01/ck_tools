function [ output_args ] = ck_rc3_visualize_network_clusterbehavior3(outdir,C2C_stat,B)
 %ck_rc3_visualize_network_clusterbehavior2
%UNTITLED Summary of this function goes here

%   Detailed explanation goes here

% C ist die C2C_stat struktur
% B .... die behavioralen DAten


assignin('base','C',C2C_stat);
assignin('base','B',B);
% fuege eine Nullzeile ein
B = add_zero(B);
show(outdir,C2C_stat,B);
 

end


function show(outdir,C2C_stat,B)
%%%%
fprintf('>>> FUNCTION START show(outdir,C2C_stat,B) \n');

fprintf('start');
handles.outdir = outdir;
handles.C2C_stat = C2C_stat;
handles.C2C = C2C_stat.C2C;
handles.B = B;
[handles.level1str, handles.level2str, handles.level3str] =  get_listboxstr_init(C2C_stat);


handles.data = C2C_stat.C2C;
handles.clusternames = handles.C2C.clusternames;

screensize = get( 0, 'Screensize' );
%screensize = [1 1 2500 1600];
screensize(4) = screensize(4)-30;
handles.current.behav_idx = 1; % der index der behavioralen Daten der gerade gzeigt wird

f = figure('Visible','off','Position',screensize,...
    'MenuBar','none',...
    'Resize','off',...
    'Toolbar','none');

%d = {'Male',52,true;'Male',40,true;'Female',25,false};
handles.axg2 = axes('Position',[0.78 0.32 0.2 0.19 ]);
handles.axg1 = axes('Position',[0.78 0.56 0.2 0.19 ]);
handles.axconn = axes('Position',[0.78 0.80 0.2 0.18 ]);
%handles.table = uitable(f,'Data',d,'Position',[0.78 0.35 0.2 0.35 ]);

%set(handles.imagesc,'ButtonDownFcn',@mydatacursorupdate);

%set(dcm_obj.CurrentDataCursor, 'Marker','o', 'MarkerFaceColor','b');

% handles.pb_next = uicontrol('Style','pushbutton',...
%               'String','next','Position',[screensize(1)-100 screensize(2)-50 40 30],...
%               'FontSize',10,...
%               'Callback',{@pb_next_Callback});
bw = 40;
bh = 25;
lbw = 100;
lbh = 130;

% listbox mit behavioralen Daten
handles.lb_behav = uicontrol('Style','listbox',...
              'String',handles.B.TG,...
              'Position',[f.Position(3)-150 20 lbw lbh],...
              'FontSize',10,...
              'Callback',{@lb_behav_Callback});

% listbox mit behavioralen Daten
handles.lb_behav_anti = uicontrol('Style','listbox',...
              'String',handles.B.TG,...
              'Position',...
              [handles.lb_behav.Position(1)-lbw-10,...
              handles.lb_behav.Position(2),lbw,lbh],...    
              'FontSize',10,...
              'min',0,'max',Inf,...
              'Callback',{@lb_behav_anti_Callback});

% ERstelle die Liste von Optionen 
strtmp = {};
if isfield(handles.C2C_stat,'C2C')
    strtmp = [strtmp; 'C2C v B'; 'C2C C'];
end
if isfield(handles.C2C_stat,'Network')
    strtmp = [strtmp; 'Net v B'; 'Net C'];
end

handles.show_var_list = {'C2C vs. B'; 'C2C C';'Net vs. B';'N C'};
handles.show_var_list = strtmp;
handles.lb_level3 = uicontrol('Style','listbox',...
              'String',handles.level3str,...
              'Position',...
              [handles.lb_behav_anti.Position(1)-lbw-10,...
              handles.lb_behav_anti.Position(2),lbw,lbh],...
              'FontSize',10,...
              'Callback',{@lb_conn_Callback});

% Listbox mit conn Daten
handles.connvarlist = get_C2C_stat_var_list_to_display(handles.outdir);

handles.lb_level2 = uicontrol('Style','listbox',...
              'String',handles.level2str,...
              'Position',...
              [handles.lb_level3.Position(1)-lbw-10,...
              handles.lb_level3.Position(2),lbw,lbh],...
              'FontSize',10,...
              'Callback',{@lb_conn_Callback});
          
handles.lb_level1 = uicontrol('Style','listbox',...
              'String',handles.level1str,...
              'Position',...
              [handles.lb_level2.Position(1)-lbw-10,...
              handles.lb_level2.Position(2),lbw,lbh],...
              'FontSize',10,...
              'Callback',{@lb_conn_Callback});
          
handles.edit_export = uicontrol('Style','edit',...
              'String','stringname',...
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

  handles.pb_diffdiff = uicontrol('Style','pushbutton',...
               'String','diag. Diff','Position',[handles.edit_export.Position(1)-10-bw*2,...
               handles.edit_export.Position(2),bw, bh],...
               'FontSize',10,...
               'TooltipString','Difference of differences',...
               'Callback',{@pb_diffdiff});
handles.edit_add_conn_to_table = uicontrol('Style','edit',...
              'String','tablename',...
              'Position',...
              [handles.pb_diffdiff.Position(1)-10-bw*2,...
               handles.pb_diffdiff.Position(2),bw*2, bh],...
              'FontSize',10,...
              'TooltipString','Table name for adding connection');
          
handles.edit_conn_new_name = uicontrol('Style','edit',...
              'String','Roi1_to_Roi2',...
              'Position',...
              [handles.edit_add_conn_to_table.Position(1)-10-bw*2,...
               handles.edit_add_conn_to_table.Position(2),bw*2, bh],...
              'FontSize',10,...
              'TooltipString','name to save in Table');
handles.pb_add_conn_to_table = uicontrol('Style','pushbutton',...
               'String','diag. Diff','Position',[handles.edit_conn_new_name.Position(1)-10-bw*2,...
               handles.edit_conn_new_name.Position(2),bw, bh],...
               'FontSize',10,...
               'TooltipString','adding the selected conn to table',...
               'Callback',{@pb_add_conn_to_table});
          

handles.edit_cm = uicontrol('Style','edit',...
              'String','[-0.5 0.5]',...
              'Position',...
              [handles.lb_level1.Position(1)-lbw-10,...
              handles.lb_level1.Position(2),bw*2,bh],...
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
          
  

         
          
% table_pos = [handles.axg1.Position(1),...
%     handles.axg1.Position(2)+handles.axg1.Position(4)+50,...
%     handles.axg1.Position(3),...
%     f.Position(4)-handles.axg1.Position(2)-handles.axg1.Position(4)-100];
%  table_pos = [0.78 0.05 0.2 0.55 ];
% handles.table = uitable('Parent', f, 'Data', cell(3,4),...
%       'Position', table_pos ,...
%       'ColumnName',{'p-value','MeanZdiff', 'Zstd', 'MeanRdiff'},...
%       'ColumnEditable',true,...
%       'ColumnWidth',{ 50 70 70 70 },...
%       'CellSelectionCallback',@(src,evnt)set(src,'UserData',evnt.Indices));
%handles.table = uitable(f,'Data',d,'Position',[0.78 0.35 0.2 0.35 ]);

%handles.table = uitable(f);
handles.table = uitable('Parent', f, 'Data', cell(9,9),...
       'Units','normalized',...
       'Position', [0.75 0.135 0.24 0.145 ],...
       'ColumnName',{'MeanRdiff','Rstd','CI_low', 'CI_high','MeanZdiff', 'Zstd', 'p-value', 't-value','df'},...
       'RowName',{'conn_g1','conn_g2','conn_dif','corr_g1','corr_g2','corr_dif','slope_g1', 'slope_g2', 'slope_dif'},...
       'ColumnEditable',true,...
    'ColumnWidth',{ 60 60 60 60 60 60 60 60 30 },...
       'CellSelectionCallback',@(src,evnt)set(src,'UserData',evnt.Indices));
% set(t,'Units','normalized');
%'ColumnWidth',{ 50 70 70 70 },...
       
% d = {'Male',52,true;'Male',40,true;'Female',25,false};
% t.Data = d;
% t.Position = [0.78 0.135 0.2 0.145 ];



%handles.G = G;
% define Axis
%assignin('base','RR12',G.RR12);
handles.axmain = axes('Position',[0.05 0.05 0.7 0.85]);
handles.imagesc = imagesc(squeeze(handles.C2C.PCC.MeanGroupDiffR(:,:)));
handles.current.name1 = 'Pearson Mean R Diff';
handles.current.name2 = 'Pearson r-value';
handles.current.imageMatrix = handles.C2C.PCC.MeanGroupDiffR;
handles.pos_old(1)=1;
handles.pos_old(2)=1;
%handles.current.imageMatrix = handles.G.RR12;

set(handles.axmain,'Units','pixels');

guidata(f,handles)

adjust_C2C_figure(handles.C2C.clusternames);
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



function adjust_C2C_figure(clusternames)
fprintf('>>> FUNCTION START adjust_C2C_figure(clusternames)\n');

handles = guidata(gcf);
areas=clusternames;
dim =length(areas);
fprintf('adjust C2C figure fkt\n');
% setze als axes die main axes
axes(handles.axmain);
set(handles.axmain,'YTick',[1:dim])
set(handles.axmain,'YTickLabel',areas)
%set(gca,'XTickLabel',{'-pi','-pi/2','0','pi/2','pi'})
set(handles.axmain,'XTick',[1:dim])
set(handles.axmain,'XTickLabel',areas')

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
fprintf('>>> FUNCTION START  myupdatefcn2) \n');

handles = guidata(gcf);
%idx = handles.current.behav_idx;
pos = get(event_obj,'Position');
% mache nur etwas wenn wir eine neue Position haben
if pos(1)~=handles.pos_old(1) || pos(2)~=handles.pos_old(2)
    
handles.pos = pos;
handles.pos_old = pos;
guidata(gcf,handles);
% loesche das alte REchteck
fprintf('before calling delete_old_rect() \n');

delete_old_rect();
fprintf('after calling delete_old_rect() \n');
if isfield(handles,'pos')
    % update fig 1
    fprintf('before call of show_subject_conn (handles.pos) in funktion  myupdatefcn2(empt,event_obj)\n');
    show_subject_conn(handles.pos);
    fprintf('nach show_subject_conn(handles.pos) in funktion  myupdatefcn2(empt,event_obj)\n');
fprintf('myupdate 3\n');
    % update fig 2
    txt = show_correlation(1,handles.pos);
    fprintf('nach show_correlation(1,handles.pos) in funktion  myupdatefcn2(empt,event_obj)\n');

    % update fig 3
    fprintf('myupdate 4\n');
    txt = show_correlation(2,handles.pos);
    fprintf('nach show_correlation(2,handles.pos) in funktion  myupdatefcn2(empt,event_obj)\n');
    fprintf('myupdate 5\n');

    % update table
    txt = show_table(handles.pos);
    fprintf('nach update table in funktion  myupdatefcn2(empt,event_obj)\n');

end
else
    fprintf('no changes in updatefcn2\n');
end
fprintf('<<< FUNCTION ENDE  myupdatefcn2) \n');

end




function txt = show_correlation(group,pos)
fprintf('>>> FUNCTION START  show_correlation(group,pos) \n');
fprintf('R1_3d uebergeben in show correlation group %.0f\n',group);
handles = guidata(gcf);
% beim wechsel von cluster zu netzwerk kann die pos zu gross sein

clusternames = handles.clusternames;
R1_3d = handles.current.R1_3d;
R2_3d = handles.current.R2_3d; 
% b = handles.B.BG1(:,handles.lb_behav.Value);
% c = squeeze(R1_3d(pos(2),pos(1),:));
% size(R1_3d)
% size(R2_3d)
% size(b)
% size(c)


if pos(1)>size(R1_3d,2) || pos(2)>size(R1_3d,1)
    pos(1) = 1;
    pos(2) = 1;
end
    

if group==1
    axes(handles.axg1);
    subjects = [1:size(R1_3d,3)];
    
    b = handles.B.BG1(:,handles.lb_behav.Value);
    assignin('base','b',b);
    c = squeeze(R1_3d(pos(2),pos(1),:));
    assignin('base','c',c);
    txt = {['Y: ',clusternames{pos(1)}],['X: ',clusternames{pos(2)}]};
    scatter(b,c,'x');
    %scatter(subjects,correll','o');
    xlabel(handles.B.TG{handles.lb_behav.Value});
    ylabel(handles.current.name2);
    
    [r,p] = corr(c,b);
    titlestr = [ clusternames{pos(2)}, clusternames{pos(1)}, '  R = ', num2str(r), ...
        '  p= ', num2str(p)];
    title(titlestr);
    handles.current.export.g1behav_idx = handles.lb_behav.Value;
    handles.current.export.g1b = b;
    handles.current.export.g1c=c;
    handles.current.export.g1titlestr = titlestr;
    handles.current.export.g1txt = txt;
    handles.current.export.g1xlabel = handles.B.TG{handles.lb_behav.Value};
    handles.current.export.g1ylabel = handles.current.name2;
    handles.current.export.filename = [clusternames{pos(1)} '_vs_' clusternames{pos(2)}];

else
    subjects = [1:size(R2_3d,3)];
    axes(handles.axg2);
    b = handles.B.BG2(:,handles.lb_behav.Value);
    c = squeeze(R2_3d(pos(1),pos(2),:));
    txt = {['Y: ',clusternames{pos(2)}],['X: ',clusternames{pos(1)}]};
    scatter(b,c,'x');
    %scatter(subjects,correll','o');
    [r,p] = corr(c,b);

    xlabel(handles.B.TG{handles.lb_behav.Value});
    ylabel(handles.current.name2);
    titlestr = [clusternames{pos(2)}, clusternames{pos(1)}, '  R = ', num2str(r), ...
        '  p= ', num2str(p)];
    title(titlestr);
    handles.current.export.g2behav_idx = handles.lb_behav.Value;
    handles.current.export.g2b = b;
    handles.current.export.g2c = c;
    handles.current.export.g2titlestr = titlestr;
    handles.current.export.g2txt = txt;
    handles.current.export.g2xlabel = handles.B.TG{handles.lb_behav.Value};
    handles.current.export.g2ylabel = handles.current.name2;
end
fprintf('show_correlation 5');
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
%tmp = size(handles.G.RR12,1);
fprintf('show_correlation 7\n');
tmp = size(R1_3d,1);
fprintf('show_correlation 8\n');

if group ==1
    handles.rectg1 = rectangle('Position',[pos(2)-0.5,pos(1)-0.5,1,1],'EdgeColor','r','LineWidth',2);
    handles.liney1g1 = rectangle('Position',[pos(1)-0.5,0,1,tmp],'EdgeColor','r','LineWidth',0.3);
    handles.liney2g1 = rectangle('Position',[pos(2)-0.5,0,1,tmp],'EdgeColor','r','LineWidth',0.3);
    handles.linex1g1 = rectangle('Position',[0,pos(2)-0.5,tmp,1.0],'EdgeColor','r','LineWidth',0.3);
    handles.linex2g1 = rectangle('Position',[0,pos(1)-0.5,tmp,1.0],'EdgeColor','r','LineWidth',0.3);
    handles.current.export.g1yp = yp;
    handles.current.export.g1linered = lr;
    handles.current.export.g1linetext = lrstr;
else
    handles.rectg2 = rectangle('Position',[pos(2)-0.5,pos(1)-0.5,1,1],'EdgeColor','r','LineWidth',2);
    handles.liney1g2 = rectangle('Position',[pos(1)-0.5,0,1,tmp],'EdgeColor','r','LineWidth',0.3);
    handles.liney2g2 = rectangle('Position',[pos(2)-0.5,0,1,tmp],'EdgeColor','r','LineWidth',0.3);
    handles.linex1g2 = rectangle('Position',[0,pos(2)-0.5,tmp,1.0],'EdgeColor','r','LineWidth',0.3);
    handles.linex2g2 = rectangle('Position',[0,pos(1)-0.5,tmp,1.0],'EdgeColor','r','LineWidth',0.3);
    handles.current.export.g2yp = yp;
    handles.current.export.g2linered = lr;
    handles.current.export.g2linetext = lrstr;
end
fprintf('show_correlation 15\n');

guidata(gcf,handles);
fprintf('<<< FUNCTION ENDE show_correlation(group,pos) \n');

end

function delete_old_rect()
fprintf('>>> FUNCTION START  delete_old_rect() -> ');
try
delete_single_old_rect('rectg1');
delete_single_old_rect('linex1g1');
delete_single_old_rect('liney1g1');
delete_single_old_rect('linex2g1');
delete_single_old_rect('liney2g1');
delete_single_old_rect('rectg2');
delete_single_old_rect('linex1g2');
delete_single_old_rect('liney1g2');
delete_single_old_rect('linex2g2');
delete_single_old_rect('liney2g2');
fprintf('... successful \n');
catch
    fprintf('... unsuccessful ... ERROR\n');
end

end

function handles = delete_single_old_rect(string)
handles = guidata(gcf);
if isfield(handles,string)
    delete(handles.(string));
end
guidata(gcf,handles)
end

function txt = show_table(pos)
fprintf('>>> FUNCTION START  5555555555555 show_table(pos)... \n');
handles = guidata(gcf);
% beim wechsel von cluster zu netzwerk kann die pos zu gross sein

clusternames = handles.clusternames;
R1_3d = handles.current.R1_3d;
R2_3d = handles.current.R2_3d; 
% b = handles.B.BG1(:,handles.lb_behav.Value);
% c = squeeze(R1_3d(pos(2),pos(1),:));
% size(R1_3d)
% size(R2_3d)
% size(b)
% size(c)
assignin('base','aa_conn_r_g1',squeeze(R1_3d(pos(2),pos(1),:)));
assignin('base','aa_conn_r_g2',squeeze(R2_3d(pos(2),pos(1),:)));
assignin('base','c1',squeeze(R1_3d(pos(2),pos(1),:)));
assignin('base','c2',squeeze(R2_3d(pos(2),pos(1),:)));
assignin('base','aa_b1',handles.B.BG1(:,handles.lb_behav.Value));
assignin('base','aa_b2',handles.B.BG2(:,handles.lb_behav.Value));
assignin('base','b1',handles.B.BG1(:,handles.lb_behav.Value));
assignin('base','b2',handles.B.BG2(:,handles.lb_behav.Value));
assignin('base','aa_BG1',handles.B.BG1);
assignin('base','aa_BG2',handles.B.BG2);

fprintf('show_table 1\n');

if pos(1)>size(R1_3d,2) || pos(2)>size(R1_3d,1)
    pos(1) = 1;
    pos(2) = 1;
end
    
fprintf('show_table 2\n');
    subjectsg1 = [1:size(R1_3d,3)];
    subjectsg2 = [1:size(R2_3d,3)];
    b1 = handles.B.BG1(:,handles.lb_behav.Value);
    c1 = squeeze(R1_3d(pos(2),pos(1),:));
    b2 = handles.B.BG2(:,handles.lb_behav.Value);
    c2 = squeeze(R2_3d(pos(1),pos(2),:));
    fprintf('show_table 3\n');
    txt = {['Y: ',clusternames{pos(2)}],['X: ',clusternames{pos(1)}]};
    
    
    

    
%     
%     % Berechnung des Unterschiedes der Korrelationen
%     [r1corr,p1corr] = corr(c1,b1);
%     fprintf('show_table 33\n');
%     [r2corr,p2corr] = corr(c2,b2);
%     fprintf('show_table 3333\n');
%     r1corr_n = 0.5*log((1+r1corr)/(1-r1corr));
%     r2corr_n = 0.5*log((1+r2corr)/(1-r2corr));
%     fprintf('show_table 333333\n');
%     SE = sqrt(1/(length(subjectsg1)-3)+1/(length(subjectsg2)-3));
%     zcorr_diff = (r1corr_n-r2corr_n)/SE;
%     fprintf('show_table 444\n');
%     p_one = normcdf(zcorr_diff);
%     fprintf('show_table 4444\n');
%     % Berechnung des unterschiedes der slopes
%     lmfg1 = fitlm(b1,c1,'linear');
%     lmfg2 = fitlm(b2,c2,'linear');
%     model_coeff_1 = table2array(lmfg1.Coefficients);
%     model_coeff_2 = table2array(lmfg2.Coefficients);
%     zlinear_diff = (model_coeff_1(2,1)-model_coeff_2(2,1))/(sqrt(model_coeff_1(2,2)^2+model_coeff_2(2,2)^2));
%     p_linear = normcdf(zlinear_diff);    
% %    xlabel(handles.B.TG{handles.lb_behav.Value});
% %    ylabel(handles.current.name2);
% fprintf('show_table 555\n');
%     table_data = handles.table.Data;
%  fprintf('show_table 55541212\n');
%     table_data{1,1} = r1corr;
% %     table_data = {'Male',52,true;'Male',40,true;'Female',25,false};
% %     fprintf('show_table 55541212\n');
% %     r1corr
%      fprintf('show_table 1111\n');
%     table_data{1,6} = p1corr;
% fprintf('show_table 555433\n');
%     table_data{2,6} = r2corr;
%     table_data{2,3} = p2corr;
% fprintf('show_table 55544\n');
%     table_data{3,2} = zcorr_diff;
%     table_data{3,3} = p_one;
%     
%     
% fprintf('show_table 3\n');
%     table_data{4,2} = zlinear_diff;
%     table_data{4,3} = p_linear;
%     fprintf('show_table 5\n');
% 
%     
% %     titlestr = [ clusternames{pos(2)}, clusternames{pos(1)}, '  R = ', num2str(r), ...
% %         '  p= ', num2str(p)];
% %     title(titlestr);
% 
% fprintf('show_table 6\n');
% yp = polyfit(b,c,1);
% lr = polyval(yp,b);
% plot(b,lr,'r-');
% lrstr = ['y = ' num2str(yp(1)) ' x + ' num2str(yp(2))];
% text(min(b),polyval(yp,min(b)),lrstr);
% hold off



% Table
%            'MeanRdiff','Rstd','CI_low', 'CI_high','MeanZdiff', 'Zstd', 'p-value', 't-value','df'},...

% conn_g1
% conn_g2
% conn_dif
% corr_g1
% corr_g2
% corr_dif
% slope_g1
% slope_g2
% slope_dif
table_data = cell(9,9);
fprintf('get_table1\n');
tmp = get_table_data_conn(c1);
assignin('base','tmp',tmp);
assignin('base','table_data',table_data);
table_data(1,1:9)= tmp;
fprintf('get_table2\n');
table_data(2,1:9) = get_table_data_conn(c2);
fprintf('get_table3\n');
table_data(3,1:9) = get_table_data_conn_dif(c1,c2);
fprintf('get_table4\n');
table_data(4,1:9) = get_table_data_corr(c1,b1);
fprintf('get_table5\n');
table_data(5,1:9) = get_table_data_corr(c2,b2);
fprintf('get_table6\n');
table_data(6,1:9) = get_table_data_corr_dif(c1,b1,c2,b2);
fprintf('get_table7\n');
table_data(7,1:9) = get_table_data_slope(c1,b1);
fprintf('get_table8\n');
table_data(8,1:9) = get_table_data_slope(c2,b2);
fprintf('get_table9x\n');
table_data(9,1:9) = get_table_data_slope_dif(c1,b1,c2,b2);
fprintf('asdfx\n');
% = ccc;
fprintf('222\n');
assignin('base','table_datas',table_data);
handles.table.Data = table_data;

 %handles.table.Data = table_data(;

fprintf('show_table 15\n');

guidata(gcf,handles);
fprintf('<<< FUNCTION ENDE show_table(pos) \n');

end


% extrahiert die Table data fuer die connecktivitiaetsanalyse fuer eine
% Gruppe
function table_data = get_table_data_conn(c1)
fprintf('start get_table_data_conn(c1)\n');
table_data = cell(1,9);
fprintf('0\n');
assignin('base','tablec1',c1);
[h,p,ci,stat] = ttest(c1);
fprintf('1\n');
table_data{1,1} = mean(c1); % mean r-value of group 1
fprintf('2\n');
table_data{1,2} = std(c1); % std of r-value of group 1
fprintf('3\n');
table_data{1,3} = ci(1); % CI of r-value of group 1
fprintf('4\n');
table_data{1,4} =  ci(2);
table_data{1,5} = mean(atanh(c1)); % z-value of group 1
table_data{1,6} = std(atanh(c1)); % std von z-value of group 1
fprintf('6\n');
%table_data{1,6} = [mean(atanh(c1))-(1.96*std(atanh(c1))/sqrt(length(c1))) mean(atanh(c1))+(1.96*std(atanh(c1))/sqrt(length(c1)))]; % CI des z-value of group 1
fprintf('7\n');
table_data{1,7} = p;% p value diff from zeros
fprintf('8\n');
table_data{1,8} = stat.tstat;% t value diff from zeros
fprintf('9\n');
table_data{1,9} = stat.df;% df
fprintf('ende get_table_data_conn(c1)\n');
end

function table_data = get_table_data_conn_dif(c1,c2)

[h,p,ci,stat] = ttest2(c1,c2);
table_data{1,1} = abs(mean(c1)-mean(c2)); % mean r-value of group diff
table_data{1,2} = stat.sd; % std of r-value of group diff
table_data{1,3} = ci(1); % CI of r-value of group 1
table_data{1,4} = ci(2);% mean(arctanh(c1)); % z-value of group 1
table_data{1,5} = 0; %std(arctanh(c1)); % std von z-value of group 1
table_data{1,6} = 0;% [mean(arctanh(c1))-(1.96*std(arctanh(c1))/sqrt(length(c1))) mean(arctanh(c1))+(1.96*std(arctanh(c1))/sqrt(length(c1)))]; % CI des z-value of group 1
table_data{1,7} = p;% p value diff from zeros
table_data{1,8} = stat.tstat;% t value diff from zeros
table_data{1,9} = stat.df;% df

end

function table_data = get_table_data_corr(c1,b1)
[r,p,ci1,ci2] = corrcoef(c1,b1);
    r1_n = 0.5*log((1+r(2,1))/(1-r(2,1)));
%    SE = sqrt(1/(length(c1)-3)+1/(length(c2)-3));

table_data{1,1} = r(2,1); % mean r-value of group 1
table_data{1,2} = 0; %std(c1); % std of r-value of group 1
table_data{1,3} = ci1(2,1); % CI of r-value of group 1
table_data{1,4} = ci2(2,1); % z-value of group 1
table_data{1,5} = r1_n; %arctanh(r); % std von z-value of group 1
table_data{1,6} = 0;% [mean(arctanh(c1))-(1.96*std(arctanh(c1))/sqrt(length(c1))) mean(arctanh(c1))+(1.96*std(arctanh(c1))/sqrt(length(c1)))]; % CI des z-value of group 1
table_data{1,7} = p(2,1);% p value diff from zeros
table_data{1,8} = 0; %stat.tstat;% t value diff from zeros
table_data{1,9} = 0;%stat.df;% df

end


function table_data = get_table_data_corr_dif(c1,b1,c2,b2)
[r1,p1,ci11,ci12] = corrcoef(c1,b1);
r1 = r1(2,1);
p1 = p1(2,1);
r1_n = 0.5*log((1+r1)/(1-r1));
[r2,p2,ci21,ci22] = corrcoef(c2,b2);
r2 = r2(2,1);
p2 = p2(2,1);

r2_n = 0.5*log((1+r2)/(1-r2));
SE = sqrt(1/(length(c1)-3)+1/(length(c2)-3));
zcorr_diff = (r1_n-r2_n)/SE;
p_one = normcdf(zcorr_diff);

table_data{1,1} = abs(r1-r2); % mean r-value of group 1
table_data{1,2} = 0; %std(c1); % std of r-value of group 1
table_data{1,3} = NaN;% [ci(1) , ci(2)]; % CI of r-value of group 1
table_data{1,4} = NaN; % ci high
table_data{1,5} = zcorr_diff; %arctanh(r); % std von z-value of group 1
table_data{1,6} = 0;% [mean(arctanh(c1))-(1.96*std(arctanh(c1))/sqrt(length(c1))) mean(arctanh(c1))+(1.96*std(arctanh(c1))/sqrt(length(c1)))]; % CI des z-value of group 1
table_data{1,7} = p_one;% p value diff from zeros
table_data{1,8} = NaN; %stat.tstat;% t value diff from zeros
table_data{1,9} = NaN;%stat.df;% df

end


    

function table_data = get_table_data_slope(c1,b1)
lmfg1 = fitlm(b1,c1,'linear');
model_coeff_1 = table2array(lmfg1.Coefficients);
%    SE = sqrt(1/(length(c1)-3)+1/(length(c2)-3));

table_data{1,1} = model_coeff_1(2,1); % mean r-value of group 1
table_data{1,2} = model_coeff_1(2,2); %std(c1); % std of r-value of group 1
table_data{1,3} = NaN; % CI of r-value of group 1
table_data{1,4} = NaN; % z-value of group 1
table_data{1,5} = NaN; %arctanh(r); % std von z-value of group 1
table_data{1,6} = NaN;% [mean(arctanh(c1))-(1.96*std(arctanh(c1))/sqrt(length(c1))) mean(arctanh(c1))+(1.96*std(arctanh(c1))/sqrt(length(c1)))]; % CI des z-value of group 1
table_data{1,7} = model_coeff_1(2,4);% p value diff from zeros
table_data{1,8} = model_coeff_1(2,3); %stat.tstat;% t value diff from zeros
table_data{1,9} = NaN;%stat.df;% df

end


function table_data = get_table_data_slope_dif(c1,b1,c2,b2)
    lmfg1 = fitlm(b1,c1,'linear');
    lmfg2 = fitlm(b2,c2,'linear');
    fprintf('1\n');
    model_coeff_1 = table2array(lmfg1.Coefficients);
    model_coeff_2 = table2array(lmfg2.Coefficients);
    zlinear_diff = (model_coeff_1(2,1)-model_coeff_2(2,1))/(sqrt(model_coeff_1(2,2)^2+model_coeff_2(2,2)^2));
    p_linear = normcdf(zlinear_diff);    
fprintf('2\n');


table_data{1,1} = abs(model_coeff_1(2,1)-model_coeff_2(2,1)); % mean r-value of group 1
table_data{1,2} = NaN; %std(c1); % std of r-value of group 1
table_data{1,3} = NaN; % CI of r-value of group 1
table_data{1,4} = NaN; % ci high
fprintf('3\n');
table_data{1,5} = zlinear_diff;% z value
table_data{1,6} = NaN;% [mean(arctanh(c1))-(1.96*std(arctanh(c1))/sqrt(length(c1))) mean(arctanh(c1))+(1.96*std(arctanh(c1))/sqrt(length(c1)))]; % CI des z-value of group 1
table_data{1,7} = p_linear;% p value diff from zeros
table_data{1,8} = NaN; %stat.tstat;% t value diff from zeros
table_data{1,9} = NaN;%stat.df;% df
fprintf('5x\n');
end
% 
% 
%     
%     
%     % Berechnung des Unterschiedes der Korrelationen
%     [r1corr,p1corr] = corr(c1,b1);
%     [r2corr,p2corr] = corr(c2,b2);
%     r1_n = 0.5*log((1+r1)/(1-r1));
%     r2corr_n = 0.5*log((1+r2corr)/(1-r2corr));
%     SE = sqrt(1/(length(subjectsg1)-3)+1/(length(subjectsg2)-3));
%     zcorr_diff = (r1corr_n-r2corr_n)/SE;
%     p_one = normcdf(zcorr_diff);
% 
%     
%     fprintf('show_table 4444\n');
%     % Berechnung des unterschiedes der slopes
%     lmfg1 = fitlm(b1,c1,'linear');
%     lmfg2 = fitlm(b2,c2,'linear');
%     model_coeff_1 = table2array(lmfg1.Coefficients);
%     model_coeff_2 = table2array(lmfg2.Coefficients);
%     zlinear_diff = (model_coeff_1(2,1)-model_coeff_2(2,1))/(sqrt(model_coeff_1(2,2)^2+model_coeff_2(2,2)^2));
%     p_linear = normcdf(zlinear_diff);    
% %    xlabel(handles.B.TG{handles.lb_behav.Value});
% %    ylabel(handles.current.name2);
% fprintf('show_table 555\n');
%     table_data = handles.table.Data;
%  fprintf('show_table 55541212\n');
%     table_data{1,1} = r1corr;
% %     table_data = {'Male',52,true;'Male',40,true;'Female',25,false};
% %     fprintf('show_table 55541212\n');
% %     r1corr
%      fprintf('show_table 1111\n');
%     table_data{1,6} = p1corr;
% fprintf('show_table 555433\n');
%     table_data{2,6} = r2corr;
%     table_data{2,3} = p2corr;
% fprintf('show_table 55544\n');
%     table_data{3,2} = zcorr_diff;
%     table_data{3,3} = p_one;
%     
%     
%     table_data{4,2} = zlinear_diff;
%     table_data{4,3} = p_linear;
%     
% end
% 
% 



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
fprintf('>>> FUNCTION START  pb_export_Callback(hObject,eventdata) \n');

handles = guidata(hObject);
% exportiert werden sollen folgende Daten
% Areaname1 Areaname2
% R1conn(:)  R2conn(:)
% Age1(:)  Age2(:)
% Table mit den Berechnungen
% abgespeichert im aktuellen Verzeichnis unter dem Namen der beiden Areale
% filename = handles.
% (varname).table = handles.table;
% (varname).Areaname1 = 
% (varname).areaname2 = 
% (varname).R1 =
% (varname).R2 = 
% (varname).Age1 =
% (varname).Age2 = ;
%exportname = handles.edit_export.String;
%assignin('base',handles.edit_export.String,handles.current.export);
%assignin('base',(varname),varname);
filename = fullfile(pwd,[handles.current.export.filename, '.mat']);
handles.current.export.table = handles.table;
export = handles.current.export;
save(filename,'export');
%assignin('base',handles.current.export.filename,handles.current.export);

end

function pb_diffdiff(hObject,eventdata)
fprintf('>>> FUNCTION START  pb_diffdiff(hObject,eventdata) \n');
handles = guidata(hObject);
M = handles.current.M;
X = M;
for i=1:size(M,1)
    for j = 1:i-1
        X(i,j)=M(i,j)-M(j,i);
    end
end
imagesc(X);
adjust_C2C_figure(handles.current.clusternames);
%assignin('base','M',M);

end


function pb_add_conn_to_table(hObject,eventdata)
fprintf('>>> FUNCTION START  pb_add_conn_to_table(hObject,eventdata) \n');
handles = guidata(hObject);
tablename = handles.edit_add_conn_to_table.String;
varname = handles.edit_conn_new_name.String;
filename = fullfile('./',[tablename + ".mat"]);
fprintf("tablename = %s   varname = %s \n",tablename, varname);
%(varname) = current.export.g1c;
g1c = handles.current.export.g1c;
g2c = handles.current.export.g2c;
Tg1 = table(g1c); %(varname);
Tg1.Properties.VariableNames{'g1c'} = varname;
Tg2 = table(g2c);
Tg2.Properties.VariableNames{'g2c'} = varname;
if exist(filename,'file') == 2
    load(filename);
    if ismember(varname, T1.Properties.VariableNames)
       T1.(varname) = []; 
       T2.(varname) = []; 
    end
    T1 = [T1 Tg1];
    T2 = [T2 Tg2];
else
    T1 = Tg1;
    T2 = Tg2;
end
save(filename,'T1','T2');
%handles.current.export.g1c
%handles.current.export.g2c

end



function pb_show_conn_Callback(hObject,eventdata)
fprintf('>>> FUNCTION START  pb_show_conn_Callback(hObject,eventdata) \n');

handles = guidata(hObject);

if ~isfield(handles,'C2C_stat')
    load(fullfile(handles.outdir,'C2C_stat.mat'));
    handles.C2C_stat = C2C_stat;

end
handles.current.showconn = 1;
handles.current.showbehav = 0;
guidata(hObject,handles);
update_all();

end
       
function pb_show_behav_Callback(hObject,eventdata)
fprintf('>>> FUNCTION START  pb_show_behav_Callback(hObject,eventdata) \n');

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
fprintf('>>> FUNCTION START  lb_behav_Callback(hObject,eventdata) \n');

update_all();
end

function lb_behav_anti_Callback(hObject,eventdata)
fprintf('>>> FUNCTION START  lb_behav_anti_Callback(hObject,eventdata) \n');

update_all();
end


function lb_conn_Callback(hObject,eventdata)
fprintf('>>> FUNCTION START  lb_conn_Callback(hObject,eventdata) \n');

update_all();

end


function pb_show_power_Callback(hObject,eventdata)
fprintf('>>> FUNCTION START show_power callback\n');
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



function show_subject_conn(pos)
fprintf('>>>>FUNCTION START: show subject conn fkt\n');
handles = guidata(gcf);
clusternames = handles.current.clusternames;
if pos(1)>length(clusternames) || pos(2)>length(clusternames)
    pos(1) = 2;
    pos(2) = 1;
end
name = [clusternames{pos(1)} ' vs. ' clusternames{pos(2)}];
fprintf('show subject conn fkt 2\n');
%'show subject conn'

c1 = 'red'; %color 1
c2 = 'blue'; %color 2
fprintf('before axes(handles.axconn) \n');

% verursachte den staendigern Aufruf von myupdate fkt in einer endlos
% schleife
axes(handles.axconn);
fprintf('after axes(handles.axconn) \n');
cla(handles.axconn);
fprintf('after cla \n');

%figure('Name',netname);
fprintf('show subject conn fkt 3\n');

R1 = handles.current.R1_3d;
R2 = handles.current.R2_3d;
Zp = handles.current.Zp; % wenn es keine Korrelationen waren ist dieser Wert nicht ausgerechnet worden
assignin('base','R1',R1);
assignin('base','R2',R2);

assignin('base','Zp',Zp);

% R1 = handles.data.R1_3d;
% R2 = handles.data.R2_3d;
% Zp = handles.data.Zp;
fprintf('show subject conn fkt 4\n');
%R1 = G.G1.R_3d;
%R2 = G.G2.R_3d;

% if ~isfield(handles,'C2C_stat_g1')
%     load(fullfile(handles.outdir,'C2C_stat_g1.mat'));
%     handles.C2C_stat_g1 = C2C_stat_g1;
% end
% if ~isfield(handles,'C2C_stat_g2')
%     load(fullfile(handles.outdir,'C2C_stat_g2.mat'));
%     handles.C2C_stat_g2 = C2C_stat_g2;
% end
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
fprintf('show subject conn fkt 5\n');

hold on
p2=barwitherr(e2,x2,y2,c2); %,'DisplayName','group1');

%legend([p1 t p2],'group1','std within Network','group2');

hold on 
xl1 = [x1(1):1:x2(end)];
yl1 = ones(length(xl1))*mean(squeeze(R1(pos(1),pos(2),:)));

plot(xl1,yl1,c1);

xl2 = [x1(1):1:x2(end)];
yl2 = ones(length(xl1))*mean(squeeze(R2(pos(1),pos(2),:)));
plot(xl2,yl2,c2);

fprintf('show subject conn fkt 6\n');
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
fprintf('show subject conn fkt 7\n');

fprintf('show subject conn fkt 7.2\n');

str = sprintf('R1 mean = %.4f  R2 mean = %.4f  R1-R2= %.4f    p = %.7f',...
    mean(R1(pos(1),pos(2),:)),...
    mean(R2(pos(1),pos(2),:)),...
    mean(R1(pos(1),pos(2),:))-mean(R2(pos(1),pos(2),:)),...
    Zp(pos(1),pos(2)));
fprintf('show subject conn fkt 8 %s\n',str);
xlabel(str);
ylabel('R-value');
str = sprintf('connectivity between areas: %s',name);
title(str,'Color','red','FontSize',12); 
assignin('base','y1',y1);
assignin('base','y2',y2);
fprintf('<<<< ENDE show subject conn fkt ende\n');
end

function varlist = get_C2C_stat_var_list_to_display(outdir)

load(fullfile(outdir,'C2C_stat.mat'));
varlist = cell(0);
if isfield(C2C_stat,'C2C')
if isfield(C2C_stat.C2C,'MeanZdiff')
    varlist{end+1} = 'MeanZdiff';
end
if isfield(C2C_stat.C2C,'MeanSyndiff')
    varlist{end+1} = 'MeanSyndiff';
end
if isfield(C2C_stat.C2C,'MeanPLVdiff')
    varlist{end+1} = 'MeanPLVdiff';
end
end
if isfield(C2C_stat,'Network')
if isfield(C2C_stat.Network,'MeanZdiff')
    varlist{end+1} = 'MeanZdiff';
end
if isfield(C2C_stat.Network,'MeanSyndiff')
    varlist{end+1} = 'MeanSyndiff';
end
if isfield(C2C_stat.Network,'MeanPLVdiff')
    varlist{end+1} = 'MeanPLVdiff';
end
end

end





function pb_visualizeCluster_Callback(hObject, eventdata)
fprintf('>>> FUNCTION START  pb_visualizeCluster_Callback(hObject, eventdata) \n');
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
   
    val = handles.pos(1);
    str = handles.Cluster{val}.name;
    
    SPM = handles.SPM_visualize.SPM;

    SPM.swd = handles.outdir ;
    VOL = handles.SPM_visualize.VOL;


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

function B = add_zero(B)
B.BG1 = [zeros(size(B.BG1,1),1),B.BG1];
B.BG2 = [zeros(size(B.BG2,1),1),B.BG2];
B.TG = ['Zero',B.TG];
end

function [M,P,R1_3d_clean,R2_3d_clean,clusternames] = estimate_main_behavioral_matrix(handles) %, R1_3d,R2_3d,idx_anti,idx, B1, B2)
fprintf('>>> FUNCTION START  estimate_main_behavioral_matrix(handles) \n');
if handles.lb_behav.Value==1 && handles.lb_behav_anti.Value==1
    
    
end
l1 = handles.lb_level1.String{handles.lb_level1.Value};
l2 = handles.lb_level2.String{handles.lb_level2.Value};
l3 = handles.lb_level3.String{handles.lb_level3.Value};

clusternames = handles.C2C_stat.(l1).clusternames;

R1_3d = handles.C2C_stat.(l1).(l2).([l3(end) '1']);
R2_3d = handles.C2C_stat.(l1).(l2).([l3(end) '2']);
R1_3d_clean = zeros(size(R1_3d));
R2_3d_clean = zeros(size(R2_3d));
% assignin('base','R1_3d',R1_3d);
%handles.lb_behav.Value;
%handles.lb_behav_anti.Value

% behavioral 1 for correlation
b1c = handles.B.BG1(:,handles.lb_behav.Value);
b1a = handles.B.BG1(:,handles.lb_behav_anti.Value);
% behavioral 2 for correlation
b2c = handles.B.BG2(:,handles.lb_behav.Value);
b2a = handles.B.BG2(:,handles.lb_behav_anti.Value);

assignin('base','b1c',b1c);
assignin('base','b1a',b1a);
assignin('base','b2c',b2c);
assignin('base','b2a',b2a);

bnames = handles.B.TG;

dim1 = size(R1_3d,1);
dim2 = size(R1_3d,2);

% assignin('base','b1c',b1c);
% assignin('base','b2c',b2c);
% assignin('base','b1a',b1a);
% assignin('base','b2a',b2a);
% 
% return

M = zeros(dim1,dim2);
P = zeros(dim1,dim2);

for i=1:dim1
    for j=1:dim2

        % entferne den Parameter of no interest (ba)
        % [b,bint,r] = regress(y,X)
        if handles.lb_behav_anti.Value>1
           % 'sdd'
            y1 = squeeze(R1_3d(i,j,:));
            y2 = squeeze(R2_3d(i,j,:));
            X1 = [ones(size(y1)) b1a];
            X2 = [ones(size(y2)) b2a];
            [b,bint,r1] = regress(y1,X1);
            [b,bint,r2] = regress(y2,X2);
            R1_3d_clean(i,j,:)=r1;
            R2_3d_clean(i,j,:)=r2;
        else
            % wenn keine Regression gerechnet werden soll
           % fprintf('nuer uebertrag\n')
            R1_3d_clean(i,j,:) = R1_3d(i,j,:);
            R2_3d_clean(i,j,:) = R2_3d(i,j,:);
        end
        if i<j
            r= squeeze(R2_3d_clean(i,j,:));
            bc = b2c;
        else
            r = squeeze(R1_3d_clean(i,j,:));
            bc = b1c;
        end
        % wenn eine Korrelation ausgerechnet und angezeigt werden soll
        if handles.lb_behav.Value>1
            if length(r)~=length(bc)
                fprintf('ungleiche Vectorlaengen der Korrelationen ...\n');
                fprintf('R2_3d_clean(%d,%d,:) laenge %d\n',i,j,length(r));
                fprintf('R1_3d_clean(%d,%d,:) laenge %d\n',i,j,length(bc));
                assignin('base','R1_3d_clean',R1_3d_clean);
                assignin('base','R2_3d_clean',R2_3d_clean);
                assignin('base','r',r);
                assignin('base','bc',bc);
            end
            [M(i,j), P(i,j)] = corr(r,bc);
        else
            %assignin('base','R1_3d_clean',R1_3d_clean);
            M(i,j) = mean(squeeze(R1_3d_clean(i,j,:)))-mean(squeeze(R2_3d_clean(i,j,:)));
            [~, P(i,j)] = ttest2(squeeze(R1_3d_clean(i,j,:)),squeeze(R2_3d_clean(i,j,:)));
        end
    end
end

end



function [ Z3d ,R_3d_clean] = get_reg_M(R,b)
fprintf('>>> FUNCTION START  get_reg_M(R,b) \n');

tic
%Z = zeros(size(R,1),size(R,2));
R_3d_clean = zeros(size(R));
Zx = zeros(size(R));

for i=1:size(R,1)
    for j=1:size(R,2)
            if i~=j
            y_noisy = squeeze(R(i,j,:));
            noise = b;
            dummy= ck_rc3_linear_regression_vector(y_noisy,noise);
            dummy(dummy<-1)=-0.99;
            dummy(dummy>1)=0.99;

            R_3d_clean(i,j,:)=dummy;
            dummy2 = ck_rc3_fisher_transformation(dummy);
            if ~isreal(dummy2)
                fprintf('not real i= %.0f   j = %.0f \n',i,j);
                error('sdf');
            end
            Z3d(i,j,:) = dummy2;
            
            end
    end
end

Z = mean(Z3d,3);
end

function [level1str, level2str, level3str] =  get_listboxstr_init(C2C_stat)
fprintf('>>> FUNCTION START  get_listboxstr_init(C2C_stat) \n');
level1str = cell(0);
level2str = cell(0);
level3str = cell(0);
f1 = fields(C2C_stat);
for i=1:length(f1)
    if isstruct(C2C_stat.(f1{i}))
        level1str{end+1} = f1{i};
    end
end

t2 = C2C_stat.(level1str{1});
f2 = fields(t2);
for i=1:length(f2)
    if isstruct(t2.(f2{i}))
        level2str{end+1} = f2{i};
    end
end

t3 = C2C_stat.(level1str{1}).(level2str{1});
f3 = fields(t3);
for i=1:length(f3)
    if isnumeric(t3.(f3{i}))
        level3str{end+1} = f3{i};
    end
end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% update Funktion
%%%% zentrale Funktion die bei den meisten clicks oder manipulationen 
%%%% aufgerufen wird und entsprechend des aktuellen Zustandes der GUI
%%%% die Ausgaben konfiguriert
function update_all(varargin)

handles = guidata(gcf);
fprintf('>>> FUNCTION START  update all\n');
% falls der string 'not main' uebergeben wird dann nicht das main fig
% ueberschreiben
if nargin==1
    str = varargin{1};
else
    str = 'all';
end

handles.current = get_main_image_matrix(handles);


if strcmp(str,'not main')
    fprintf('main window will not be updated');
else
    %axes(handles.axmain);
    imagesc(handles.current.M);
    adjust_C2C_figure(handles.current.clusternames);
end
guidata(gcf,handles);


if isfield(handles,'pos')
    % update fig 1
    show_subject_conn(handles.pos);
    
    % update fig 2
    txt = show_correlation(1,handles.pos);
    
    % update fig 3
    txt = show_correlation(2,handles.pos);
end

end


function current = get_main_image_matrix(handles)
fprintf('>>> FUNCTION START  get main image matrix\n');
%handles = guidata(gcf);
if ~isfield(handles,'C2C_stat')
    load(fullfile(handles.outdir,'C2C_stat.mat'));
    handles.C2C_stat = C2C_stat;
end

handles.current.name1 = handles.level3str; % m
handles.current.name2 = 'Platzhalter'; % m
if strcmp(handles.level3str(end),'R')
    handles.current.name2 = 'Pearson r-value';
end
            %[M,Zp,R1_3d_clean,R2_3d_clean] = estimate_main_connectivity_matrix(handles);
            [M,Zp,R1_3d_clean,R2_3d_clean,clusternames] = estimate_main_behavioral_matrix(handles);

            handles.current.M = M;
            handles.current.Zp = Zp;
            handles.current.R1_3d = R1_3d_clean;
            handles.current.R2_3d = R2_3d_clean;
            handles.current.clusternames = clusternames;
            
    

%guidata(gcf,handles);
current = handles.current;
assignin('base','current_at_end_of_get_main_image_matrix',current);
end


