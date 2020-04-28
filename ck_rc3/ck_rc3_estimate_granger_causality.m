function cGCt = ck_rc3_estimate_granger_causality(outdir,datadir,filenames,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Es stehen 2 Methoden zur Berechnung der Granger Causality zur Verfügung
% 1. Berechung der conditional Granger causality in the time domain
% 2. Berechnung der conditional Granger causality in the frequency domain
%
% Wichtige paper:
% Geweke 1984 Measures of conditional linear depedence and feedback between
% time series
% --- hier werden die Mathematischen Grundlagen beschrieben

% Liao 2010 Evaluating the effective connectivity of resting state networks
% using conditional ggranger causality
% --- Anwendung der conditional Granger causality auf resting state
% --- Networks ohne wesentliche neue methodische Aspekte

% Chen 2006 Frequency decomposition of conditional Granger causality and
% application to multivariate neural field potential data
% --- Erweiterung der CGC von der zeit in die Frequenz

% Zhou 2009 Analyzing Brain Networks with PCA and Conditional Granger
% Causality
% --- Anwendung der CGC in der frequency domain auf echte und
% --- simulierte Daten. Die theorie ist besser bei Chen nachzulesen
% --- aber die Autoren bieten hier Formeln für simulierte Daten
% --- an (welche im Program ck_rc_simulate_CGCA_data implementiert sind)
% --- Hiermit können die Ergebnisse die dargestellt werden mit den
% --- eigenen Methoden zum Zwecke der Qualitätssicherung verglichen werden.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if options.cancel
    return
end
load(fullfile(outdir,'parameter.mat'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Berechne die Statistik von einer Gruppe mittels
    % der vorimplementierten mtrial methode
    % aus der Granger causality toolbox
    if options.onegroup
        fprintf('for one group...\n');
        freq1=options.freq(1);
        step=options.prep.freqiteration;
        freq2=options.freq(2);
        [Vm,Pm] = ck_rc_multisubject_selection();
        
        for i=freq1:step:freq2
            options.freq(1)=i;
            options.freq(2)=i+step;
            X=configure_one_group_matrix(Vm,Pm,outdir,options);
            assignin('base','Xonegorup',X);
            %        h = waitbar(0,'please wait');
            %for i=1:size(Pm,1)
            %   waitbar(i/size(Pm,1));
            [bic,aic] = cca_find_model_order_mtrial(X,size(X,1),size(Vm{1,1},1),1,15);
            fprintf('estimated model order_: %.f / %.f\n',bic,aic);
            cGCt.ret=cca_granger_regress_mtrial(X,size(Pm,1),size(Vm{1,1},1),options.nlagst,1);
            cGCt.nodenames=options.Cluster_sel_name;
            assignin('base','cGCt_mtrial',cGCt);
            save(fullfile(outdir,['cGCt_groupA_mtrial' num2str(options.freq(1)) '.mat']),'cGCt');
        end
        %cGCt{i,1}.nodenames=options.Cluster_sel_name;
        %end
        %       close(h);
    end
    % Berechne die Statistik von einer Gruppe mittels
    % über eine separate Statisik jedes einzelnen Probanden
    if options.onegroupsingle
        fprintf('for one group single...\n');
        [Vm1,Pm1] = ck_rc_multisubject_selection();
        cGCt_groupA=estimate_one_group_granger_causality(Vm1,Pm1,outdir,options);
        assignin('base','cGCt_groupA',cGCt_groupA);
        save(fullfile(outdir,'cGCt_groupA.mat'),'cGCt_groupA');
        %         [Vm,Pm] = ck_rc_multisubject_selection();
        %         h = waitbar(0,'please wait');
        %         for i=1:size(Pm,1)
        %             waitbar(i/size(Pm,1));
        %             % lade die Daten
        %             [Data.Y, Data.XYZ] = spm_read_vols(Vm{i,1});
        %             %speichern der Daten
        %             filenamep = ['Data_' Pm{i,1} '.mat'];
        %             save(fullfile(outdir,filenamep),'Data');
        %             % Fuelle den Cluster und lade Ihn anschließend
        %             ck_rc_fill_cluster(outdir,filenamep);
        %             load(fullfile(outdir,'Cluster.mat'));
        %             %           load(fullfile(outdir,filenamep));
        %             %Cluster
        %             Cluster=Cluster(options.Cluster_sel_num,1);
        %             X = preprocessing(Data,Cluster,Vm{i,1}(1,1),outdir,filenamep,options,Pm{i,1});
        %             if i==1
        %                 [bic,aic] = cca_find_model_order(X,1,15);
        %                 fprintf('estimated model order: %.f / %.f\n',bic,aic);
        %             end
        %             cGCt{i,1}=estimate(X,outdir,Pm{i,1},options);
        %             cGCt{i,1}.nodenames=options.Cluster_sel_name;
        %         end
        %         close(h);
    end
    if options.twogroups
        fprintf('for two group...\n');
        [Vm1,Pm1] = ck_rc_multisubject_selection();
        [Vm2,Pm2] = ck_rc_multisubject_selection();
        cGCt_groupA=estimate_one_group_granger_causality(Vm1,Pm1,outdir,options);
        cGCt_groupB=estimate_one_group_granger_causality(Vm2,Pm2,outdir,options);
        assignin('base','cGCt_groupA',cGCt_groupA);
        assignin('base','cGCt_groupB',cGCt_groupB);
    end
    if options.pairedgroups
        fprintf('for paired groups...\n');
        [Vm1,Pm1] = ck_rc_multisubject_selection();
        [Vm2,Pm2] = ck_rc_multisubject_selection();
        assignin('base','Vm2_tmp',Vm2);
        assignin('base','Vm1_tmp',Vm1);
        cGCt_groupA=estimate_one_group_granger_causality(Vm1,Pm1,outdir,options);
        cGCt_groupB=estimate_one_group_granger_causality(Vm2,Pm2,outdir,options);
        assignin('base','cGCt_groupA',cGCt_groupA);
        assignin('base','cGCt_groupB',cGCt_groupB);
        cGCt_stat=ck_rc_estimate_cGCt_paired_group_statistik(cGCt_groupA, cGCt_groupB,options.pt);
        assignin('base','cGCt_stat',cGCt_stat);
        
    end



%assignin('base','cGCt',cGCt);
%assignin('base','X',X);
%save(fullfile(outdir,'cGCt.mat'),'cGCt');
%ck_rc_analyse_cGCt(cGCt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function cGCt=estimate_one_group_granger_causality(outdir,datadir,filenames,options)
%changed 20.05.2013
Data=0;
filenamep='x';
h = waitbar(0,'please wait');
%ck_rc_fill_all_clusters(outdir,Vm,Pm,'forced');

for i=1:size(filenames,1)
    waitbar(i/length(filenames));
    % lade die Daten
    %     [Data.Y, Data.XYZ] = spm_read_vols(Vm{i,1}); %ck 20.05.2013
    %speichern der Daten
    %filenamep = ['Data_' Pm{i,1} '.mat'];
    %save(fullfile(outdir,filenamep),'Data');
    % Fuelle den Cluster und lade Ihn anschließend
    %ck_rc_fill_cluster(outdir,filenamep);
    load(fullfile(datadir,[filenames{i} ]));
    %           load(fullfile(outdir,filenamep));
    %Cluster
    Cluster=Cluster(options.Cluster_sel_num,1);
    
    X = preprocessing(Cluster,outdir,filenamep,options,Pm{i,1});
    
    if i==1
        
        [bic,aic] = cca_find_model_order(X,1,15);
        
        fprintf('estimated model order: %.f / %.f\n',bic,aic);
    end
    cGCt{i,1}=estimate(X,outdir,Pm{i,1},options);
    cGCt{i,1}.nodenames=options.Cluster_sel_name;
end
close(h);



function X = preprocessing(Cluster,outdir,filenamep,options,znames)
% extract time series
X=zeros(size(Cluster{1,1}.Y_mean,1),size(Cluster,1));

for i=1:size(Cluster,1)
    X(:,i)=Cluster{i,1}.Y_mean;
end
%X=X'; % changed 03.06.2013 Spaltenvector sollte es sein
% detrending
if options.prep.detrend
    for i=1:size(X,1)
        X(i,:)= detrend(X(i,:));
    end
end
% removemean
if options.prep.removemean
    [X,m,e] = cca_rm_temporalmean(X,0);
end
% removemeanstd
if options.prep.removemeanstd
    [X,m,e] = cca_rm_temporalmean(X,1);
end
% differencing
if options.prep.differencing
%    [X] = cca_diff(X);    
    [X] = cca_diff(X');    %changed by sb 11.4.14
    
end

if options.prep.filter
    % freqeny band pass filter
    % options.freq
    assignin('base','X',X);
    TR=options.TR;
    Fpass1=options.freq(1);
    Fstop1=Fpass1-0.005;
    Fpass2=options.freq(2);
    Fstop2=Fpass2+0.005;
    for i=1:size(X,1)
        X(i,:)=ck_rc_bandpass_filter_slow(X(i,:),Fstop1,Fpass1,Fpass2,Fstop2,TR);
    end
    assignin('base','Xfilter',X);
end


function cGCt = estimate(X,outdir,znames,options)

%%%%%%%%%%%%%%%% CGC in time domain %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options.EstimateCGCt
    cGCt.ret=ck_rc_cca_granger_regress(X,options.nlagst);
    cGCt.name=znames;
    [cGCt.PR,cGCt.q] = cca_findsignificance(cGCt.ret,options.pt,options.cort);
    cGCt.ret.cdensity=cca_causaldensity(cGCt.ret.gc);
    cGCt.ret.cflow=cca_causalflow(cGCt.ret.gc);
end
if options.EstimateGCautonomy
    cGCt.aut.ret = cca_autonomy_regress(X,options.nlagst);
    cGCt.aut.name=znames;
    [cGCt.aut.PR,cGCt.aut.q] = cca_findsignificance_autonomy(cGCt.aut.ret,options.pt,options.cort);
end
if options.EstimateCGCf
    %cGCt.ret = estimate_CGCtX,outdir,options,znames);
    %[GW,COH,pp,waut,cons]=cca_pwcausal(X,Nr,Nl,nlags,Fs,freq,STATFLAG)
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% CGC in time domain %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function ret = estimate_CGCt(X,outdir,options,zname)
% load(fullfile(outdir,'parameter.mat'));
% %V_org = V;
% %idx=0;
% nlags=options.nlagst;
% %%%%%%%%%%%%%%%%%%%%%%%
% % erstelle eine mxn Matrix X der Clusterdaten
% % mit m... Anzahl der Messpunkte
% %     n... Anzahl der Cluster
% % M=zeros(size(Cluster{1,1}.Y_mean,1),size(Cluster,1));
% % for i=1:size(Cluster,1)
% %     M(:,i)=Cluster{i,1}.Y_mean;
% % end
% % M=M';
% % assignin('base','M',M);
% % rufe nun die Toolbox zur Berechnung der CGC auf
% [ret] = cca_granger_regress(X,nlags);
% %[PR,q] = cca_findsignificanc(ret,options.pt,options.cort);
%
% assignin('base','ret',ret);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% CGC in frequency domain %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function estimate_CGCf(Data,Cluster,V,outdir,datafilename,options,zname)
% load(fullfile(outdir,'parameter.mat'));
% fprintf('Berechnung der conditianl Granger causality in frequency domain\n');
% idx = 0;

function X= configure_one_group_matrix(Vm,Pm,outdir,options)


Data=0;
filenamep='x';
h = waitbar(0,'please wait');
ck_rc_fill_all_clusters(outdir,Vm,Pm);

for i=1:size(Pm,1)
    waitbar(i/size(Pm,1));
    
    load(fullfile(outdir,['Cluster_' Pm{i,1} '.mat']));
    Cluster=Cluster(options.Cluster_sel_num,1);
    if i==1
        X = preprocessing(Data,Cluster,Vm{i,1}(1,1),outdir,filenamep,options,Pm{i,1});
        [bic,aic] = cca_find_model_order(X,1,15);
        fprintf('estimated model order: %.f / %.f\n',bic,aic);
    else
        Y = preprocessing(Data,Cluster,Vm{i,1}(1,1),outdir,filenamep,options,Pm{i,1});
        X = [X Y];
    end
end
close(h);





