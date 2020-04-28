function ck_rc3_estimate_connectivity(clusterdir,outdir,filelist,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Es stehen 2 Methoden zur Berechnung der Konnektivität zur Verfügung
% 1. Berechung von Correlationen
% 2. Berechnung von Coherencen
%
% Die Berechnung der Connectivität nach Correlationen ist nachgebildet 
% nach dem paper :
%     "Whole brain functional connectivity in the early blind"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation of inter-regional Pearson's correlations
% Regional mean time series were estimated by averaging the time 
% series of all voxels in this region 
% (Salvador et al., 2005aGo, bGo; Achard et al., 2006Go; Liang et al., 2006Go). 
% The Pearson's correlation coefficients were computed between each pair 
% of brain regions for each subject. For further statistical analysis, 
% a Fisher's r-to-z transformation z = 0.5 x log[(1 + r)/(1 – r)] was 
% applied to improve the normality of the correlation coefficients. 
% The individual z scores were entered into a one-sample two-tailed t-test 
% to determine if the two brain regions show significant functional 
% connectivity within each group. They were also entered into a two-sample 
% two-tailed t-test to determine if the functional connectivities were 
% significantly different between the two groups.
% 
% A t-test was performed for all the 6670 (116 x 115/2) 
% functional connectivities, so a correction for multiple comparisons 
% was strictly necessary. The false discovery rate (FDR) approach was 
% applied to find a threshold that would restrict the expected proportion 
% of type I errors to lower than 0.05 (Benjamini and Yekutieli, 2001Go; 
% Salvador et al., 2005aGo). In this study, we identified the significant 
% differences in functional connectivities between the blind and sighted 
% subjects according to the following two criteria: (a) the z values were 
% significantly different from zero at least in one group at P < 0.05 
% (one-sample two-tailed t-test; FDR corrected); (b) the z scores were 
% significantly different between the two groups at P < 0.05 (two-sample 
% two-tailed t-test; FDR corrected).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Die Berechnung der Connectivität nach Coherencen ist nachgebildet 
% nach dem paper :
%     "Measuring interregional functional connectivity using 
%      coherence and partial coherence analyses of fMRI data"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%options = ck_rc2_estimate_connectivity_options();
%assignin('base','options',options);
%assignin('base','filelist',filelist);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fprintf('performing MultiSubjectAnalysis...\n');

% Schleife über die ausgewählten Probanden
h = waitbar(0,'please wait');

if options.Cluster2Cluster
    load(fullfile(clusterdir,'Cluster.mat'));
    ClusterZ=zeros(size(Cluster,1),size(Cluster,1),1);
    save(fullfile(outdir,'ClusterZ.mat'),'ClusterZ');
end
if options.Voxel2Voxel
    %hole zusatzinfos wie z.B. die c1 threshold oder
    % auch die Multiprozessoranzahl
end
for i=1:length(filelist)
    
    waitbar(i/length(filelist));
    % lade die Daten
    %[Data.Y, Data.XYZ] = spm_read_vols(Vm{i,1});
    %ck_rc_fill_cluster(outdir,Data);
    
    %load(fullfile(outdir,'Cluster.mat'));
    datafile = fullfile(outdir,filelist{i});
    postfix = get_postfix(datafile);
    % lade die Cluster Struktur
    load(datafile);
    %Data
    %Berechne die Korrelationen
    if options.EstimateCorrelations
        if options.useY
             Cluster = estimate_corr(Cluster,outdir,options,postfix,'');
        end
        if options.useYc1
            Cluster = estimate_corr(Cluster,outdir,options,postfix,'c1');
        end
        if options.useYc2
            Cluster = estimate_corr(Cluster,outdir,options,postfix,'c2');
        end
    end
    if options.EstimateCoherence
        if options.useY
            Cluster = estimate_cohe(Cluster,outdir,options,'');
        end
        if options.useYc1
            Cluster = estimate_cohe(Cluster,outdir,options,'c1');
        end
        if options.useYc2
            Cluster = estimate_cohe(Cluster,outdir,options,'c2');
        end
        %Cluster = estimate_cohe(Data,Cluster,Vm{i,1}(1,1),outdir,options,filelist{i});
        
    end
    if options.EstimatePartialCoherence
        Cluster = estimate_part_cohe(Data,Cluster,Vm{i,1}(1,1),outdir,options,filelist{i});
    end
    if options.EstimateTransinformation
        %Cluster = estimate_transinformation(Data,Cluster,Vm{i,1}(1,1),outdir,options,filelist{i});
        if options.useY
            Cluster = estimate_transinformation(Cluster,outdir,options,postfix,'');
        end
        if options.useYc1
            Cluster = estimate_transinformation(Cluster,outdir,options,postfix,'c1');
        end
        if options.useYc2
            Cluster = estimate_transinformation(Cluster,outdir,options,postfix,'c2');
        end

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SYNCHRONIZATION
    if options.EstimateSynchronization
        %Cluster = estimate_synchronization(Data,Cluster,Vm{i,1}(1,1),outdir,options,filelist{i});
        if options.useY
            Cluster = estimate_synchronization(Cluster,outdir,options,postfix,'');
        end
        if options.useYc1
            Cluster = estimate_synchronization(Cluster,outdir,options,postfix,'c1');
        end
        if options.useYc2
            
            Cluster = estimate_synchronization(Cluster,outdir,options,postfix,'c2');
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COMPLEXITY
    if options.EstimateComplexity
        %Cluster = estimate_synchronization(Data,Cluster,Vm{i,1}(1,1),outdir,options,filelist{i});
        if options.useY
            Cluster = estimate_complexity(Cluster,outdir,options,postfix,'');
        end
        if options.useYc1
        Cluster = estimate_complexity(Cluster,outdir,options,postfix,'c1');
        end
        if options.useYc2
            
        Cluster = estimate_complexity(Cluster,outdir,options,postfix,'c2');
        end
    end
    save(fullfile(outdir,['Cluster_' postfix '.mat']),'Cluster','-v7.3');
end
close(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Cluster = estimate_corr(Cluster,outdir,options,zname,sel)
%load(fullfile(outdir,'parameter.mat'));
fprintf('Berechnung der Correlationen\n');

r_string = ['R' sel];
p_string = ['P' sel];
z_string = ['Z' sel];
y_string = ['Y' sel '_mean'];
idx=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options.Cluster2Cluster
    % warning off
    if size(Cluster,1)>1
        if isfield(Cluster{1,1}, y_string) %ADDED by sb 150518
            mean_tc = zeros(length(Cluster),length(Cluster{1,1}.(y_string)));
            for j=1:size(Cluster,1)
                mean_tc(j,:) = Cluster{j,1}.(y_string);
            end
            for j=1:length(Cluster)
                idx=idx+1;
                [R{idx},P{idx}] = ck_rc3_get_pearsons_correlations(Cluster{j,1}.(y_string),mean_tc);
                Z{idx} = ck_rc3_fisher_transformation(R{idx});
                % Speichere die Analyse in der Clusterstruktur
                Cluster{j,1}.Cluster2Cluster.(r_string) = R{idx};
                Cluster{j,1}.Cluster2Cluster.(p_string) = P{idx};
                Cluster{j,1}.Cluster2Cluster.(z_string) = Z{idx};
            end
        end
    else
        msgbox('you must define at least 2 Clusters for Cluster2Cluster analysis');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%
if options.Cluster2Voxel 
    for j=1:size(Cluster,1)
        idx=idx+1;
        [R{idx},P{idx}] = ck_rc3_get_pearsons_correlations(Cluster{j,1}.Y_mean,Data.Y);
        Z{idx} = ck_fisher_transformation(R{idx});
        % Speichere die Analyse in der Clusterstruktur
        Cluster{j,1}.Cluster2Voxel.R = R{idx};
        Cluster{j,1}.Cluster2Voxel.P = P{idx};
        Cluster{j,1}.Cluster2Voxel.Z = Z{idx};
       
        % speichere Z als Bild
        V.fname = fullfile(outdir,[ 'Z_map_C2V_' zname '_' Cluster{j,1}.name '_.nii']);
        V.descrip = ['connectivity Z-map Cluster2Voxel ' Cluster{j,1}.name];
        rmfield(V,'pinfo');
        % setze die Genauigkeit der Grauwerte
        V.dt= [16 0];
        spm_write_vol(V,Z{idx});
        
        %%%% smoothing des Bildes wenn gewünscht
        if options.C2V_smoothing
            % smoothe das Z-Bild
            [p f e] = fileparts (V.fname);
            sfilename = fullfile(p,['s' f e]);
            spm_smooth(V,sfilename,options.C2V_kernel);
            V.fname = sfilename;
        end
        % Begrenzung des Bildes auf die graue Substanz
        if options.C2V_c1mask
            [p f e] = fileparts (V.fname);
            tfilename = fullfile(p,['t' f e]);
            % einlesen des individuellen c1 bildes
            [p f e] = fileparts(V_org.fname);
            % suche c1 Bild
            c1 = dir(fullfile(p,'c1*.nii'));
            if size(c1,1)~=1
                errordlg('error during searching c1 file for masking');
                return
            else
                c1_filename = fullfile(p,c1(1,1).name);
                V_c1 = spm_vol(c1_filename);
                [Y_c1, XYZ_c1] = spm_read_vols(V_c1);
                dummy=Y_c1>options.C2V_threshold;
                [Z{idx}, XYZ_Z] = spm_read_vols(V);
                Z_new = Z{idx}.*dummy;
                V.fname = tfilename;
                %assignin('base','Z_org',Z{idx});
                %assignin('base','Z_new',Z_new);
                %assignin('base','V',V);
                %assignin('base','Y_c1',Y_c1);
                spm_write_vol(V,Z_new);
            end
        end
        
        
    end
end


if options.Voxel2Voxel
    %ck_rc2_get_pearsons_correlations();
    %disp('Voxel2Voxel is not implemented jet');
    voxel_num = size(Data.Y,1)*size(Data.Y,2)*size(Data.Y,3);
    
    % einlesen des individuellen c1 bildes
    [p f e] = fileparts(V_org.fname);
    % suche c1 Bild
    c1 = dir(fullfile(p,'c1*.nii'));
    if size(c1,1)~=1
        errordlg('error during searching c1 file for masking');
        return
    else
        c1_filename = fullfile(p,c1(1,1).name);
        V_c1 = spm_vol(c1_filename);
        [Y_c1, XYZ_c1] = spm_read_vols(V_c1);
        Y_c1_bool=Y_c1>options.C2V_threshold;
    end
    R=zeros(size(Data.Y,1),size(Data.Y,2),size(Data.Y,3));
    Rp=zeros(size(Data.Y,1),size(Data.Y,2),size(Data.Y,3));
    Rn=zeros(size(Data.Y,1),size(Data.Y,2),size(Data.Y,3));
    %P=zeros(size(Data.Y,1),size(Data.Y,2),size(Data.Y,3));
    Z=zeros(size(Data.Y,1),size(Data.Y,2),size(Data.Y,3));
    Zp=zeros(size(Data.Y,1),size(Data.Y,2),size(Data.Y,3));
    Zn=zeros(size(Data.Y,1),size(Data.Y,2),size(Data.Y,3));
    Zstd=zeros(size(Data.Y,1),size(Data.Y,2),size(Data.Y,3));
    Rstd=zeros(size(Data.Y,1),size(Data.Y,2),size(Data.Y,3));
    num_c1_voxel=sum(sum(sum(Y_c1_bool)));
    %fprintf('vor der Hauptschleife\n');
    

    % Schleife ?ber alle Voxel die grauer Substanz entsprechen
    %hh=waitbar(0,'Voxel2Voxel wait...');
    for j=1:size(Data.Y,1)
        for k=1:size(Data.Y,2)
            %waitbar(k/size(Data.Y,2));
	     fprintf('j= %.0f    k=%.0f\n',j,k);
            for l=1:size(Data.Y,3)
                if Y_c1_bool(j,k,l)
                    %fprintf('Main\n');
                    % nun wird von diesem Voxel die fc zu jedem anderen
                    % Voxel in der grauen Substanz des Gehirns berechnet
                    % daraus wird eine Karte erstellt
                    % und der Mittelwert aller Voxel wird zum neuen Wert
                    % die stderr wird in einer eigenen Karte gespeichert
                    tc=squeeze(Data.Y(j,k,l,:));
                    [Rt,~] = ck_rc3_get_pearsons_correlations(tc,Data.Y);
                    %assignin('base','Rt',Rt);
                    %assignin('base','Pt',Pt);
                    %save('H:\data\Vesti\auswertung\test\Rt.mat','Rt');
%                    assignin('base','Pt',Pt);
                    Zt = ck_fisher_transformation(Rt);
                    Rt=Rt.*Y_c1_bool;
                    Zt=Zt.*Y_c1_bool;
                    Y_p_bool=Zt>0;
		            num_p_voxel=sum(sum(sum(Y_p_bool)));
                    Y_n_bool=Zt<0;
		            num_n_voxel=sum(sum(sum(Y_n_bool)));
                    R(j,k,l)=sum(sum(sum(Rt)))*100/num_c1_voxel;
                    Rp(j,k,l)=sum(sum(sum(Rt(Rt>0))))*100/(num_p_voxel);
                    Rn(j,k,l)=sum(sum(sum(Rt(Rt<0))))*100/(num_n_voxel);
                    S=reshape(Rt(Y_c1_bool),sum(sum(sum(Y_c1_bool))),1);
                    Rstd(j,k,l)=std(S);
                    Z(j,k,l)=sum(sum(sum(Zt)))*100/num_c1_voxel;
                    Zp(j,k,l)=sum(sum(sum(Zt(Zt>0))))*100/(num_p_voxel);
                    Zn(j,k,l)=sum(sum(sum(Zt(Zt<0))))*100/(num_n_voxel);
                    S=reshape(Zt(Y_c1_bool),sum(sum(sum(Y_c1_bool))),1);
                    Zstd(j,k,l)=std(S);
                    
                    
                end
            end
            
        end
        % zwischenspeicherung
        % speichere Z als Bild
        V.fname = fullfile(outdir,[ 'Z_map_V2V_' zname '.nii']);
        V.descrip = ['connectivity Z-map Voxel2Voxel ' zname];
        rmfield(V,'pinfo');
        % setze die Genauigkeit der Grauwerte
        V.dt= [16 0];
        spm_write_vol(V,Z);
        V.dt= [16 0];
        spm_write_vol(V,Z);
        V.fname = fullfile(outdir,[ 'Zp_map_V2V_' zname '.nii']);
        V.descrip = ['connectivity Zp-map Voxel2Voxel ' zname];
        rmfield(V,'pinfo');
        % setze die Genauigkeit der Grauwerte
        V.dt= [16 0];
        spm_write_vol(V,Zp);
        V.fname = fullfile(outdir,[ 'Zn_map_V2V_' zname '.nii']);
        V.descrip = ['connectivity Zn-map Voxel2Voxel ' zname];
        rmfield(V,'pinfo');
        % setze die Genauigkeit der Grauwerte
        V.dt= [16 0];
        spm_write_vol(V,Zn);
        % speichere R als Bild
        V.fname = fullfile(outdir,[ 'R_map_V2V_' zname '.nii']);
        V.descrip = ['connectivity R-map Voxel2Voxel ' zname];
        rmfield(V,'pinfo');
        % setze die Genauigkeit der Grauwerte
        V.dt= [16 0];
        spm_write_vol(V,R);
        V.fname = fullfile(outdir,[ 'Rp_map_V2V_' zname '.nii']);
        V.descrip = ['connectivity Rp-map Voxel2Voxel ' zname];
        rmfield(V,'pinfo');
        % setze die Genauigkeit der Grauwerte
        V.dt= [16 0];
        spm_write_vol(V,Rp);
        V.fname = fullfile(outdir,[ 'Rn_map_V2V_' zname '.nii']);
        V.descrip = ['connectivity Rn-map Voxel2Voxel ' zname];
        rmfield(V,'pinfo');
        % setze die Genauigkeit der Grauwerte
        V.dt= [16 0];
        spm_write_vol(V,Rn);
        
    end
    close(hh);
    % speichere Z als Bild
    V.fname = fullfile(outdir,[ 'Z_map_V2V_' zname '.nii']);
    V.descrip = ['connectivity Z-map Voxel2Voxel ' zname];
    rmfield(V,'pinfo');
    % setze die Genauigkeit der Grauwerte
    V.dt= [16 0];
    spm_write_vol(V,Z);
    V.fname = fullfile(outdir,[ 'Zp_map_V2V_' zname '.nii']);
    V.descrip = ['connectivity Zp-map Voxel2Voxel ' zname];
    rmfield(V,'pinfo');
    % setze die Genauigkeit der Grauwerte
    V.dt= [16 0];
    spm_write_vol(V,Zp);
    V.fname = fullfile(outdir,[ 'Zn_map_V2V_' zname '.nii']);
    V.descrip = ['connectivity Zn-map Voxel2Voxel ' zname];
    rmfield(V,'pinfo');
    % setze die Genauigkeit der Grauwerte
    V.dt= [16 0];
    spm_write_vol(V,Zn);
    % speichere Zstd als Bild
    V.fname = fullfile(outdir,[ 'Zstd_map_V2V_' zname '.nii']);
    V.descrip = ['connectivity Zstd-map Voxel2Voxel ' zname];
    rmfield(V,'pinfo');
    % setze die Genauigkeit der Grauwerte
    V.dt= [16 0];
    spm_write_vol(V,Zstd);
    % speichere R als Bild
    V.fname = fullfile(outdir,[ 'R_map_V2V_' zname '.nii']);
    V.descrip = ['connectivity R-map Voxel2Voxel ' zname];
    rmfield(V,'pinfo');
    % setze die Genauigkeit der Grauwerte
    V.dt= [16 0];
    spm_write_vol(V,R);
    V.fname = fullfile(outdir,[ 'Rp_map_V2V_' zname '.nii']);
    V.descrip = ['connectivity Rp-map Voxel2Voxel ' zname];
    rmfield(V,'pinfo');
    % setze die Genauigkeit der Grauwerte
    V.dt= [16 0];
    spm_write_vol(V,Rp);
    V.fname = fullfile(outdir,[ 'Rn_map_V2V_' zname '.nii']);
    V.descrip = ['connectivity Rn-map Voxel2Voxel ' zname];
    rmfield(V,'pinfo');
    % setze die Genauigkeit der Grauwerte
    V.dt= [16 0];
    spm_write_vol(V,Rn);
    % speichere Rstd als Bild
    V.fname = fullfile(outdir,[ 'R_std_map_V2V_' zname '.nii']);
    V.descrip = ['connectivity Rstd-map Voxel2Voxel ' zname];
    rmfield(V,'pinfo');
    % setze die Genauigkeit der Grauwerte
    V.dt= [16 0];
    spm_write_vol(V,Rstd);
    
    Zf = ck_fisher_transformation(R);
    Zfp = ck_fisher_transformation(Rp);
    Zfn = ck_fisher_transformation(Rn);
        % speichere Z als Bild
    V.fname = fullfile(outdir,[ 'Zf_map_V2V_' zname '.nii']);
    V.descrip = ['connectivity Zf-map Voxel2Voxel ' zname];
    rmfield(V,'pinfo');
    % setze die Genauigkeit der Grauwerte
    V.dt= [16 0];
    spm_write_vol(V,Zf);
    V.fname = fullfile(outdir,[ 'Zfp_map_V2V_' zname '.nii']);
    V.descrip = ['connectivity Zfp-map Voxel2Voxel ' zname];
    rmfield(V,'pinfo');
    % setze die Genauigkeit der Grauwerte
    V.dt= [16 0];
    spm_write_vol(V,Zfp);
    V.fname = fullfile(outdir,[ 'Zfn_map_V2V_' zname '.nii']);
    V.descrip = ['connectivity Zfn-map Voxel2Voxel ' zname];
    rmfield(V,'pinfo');
    % setze die Genauigkeit der Grauwerte
    V.dt= [16 0];
    spm_write_vol(V,Zfn);

end

%assignin('base','R',R{idx});
%assignin('base','P',P{idx});
%assignin('base','Z',Z{idx});
%save(fullfile(outdir,['Cluster_' zname '.mat']),'Cluster','-v7.3');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% COHERENCE ANALYSIS %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nach folgendem Paper
% Measuring interregional functional connectivity using coherence and
%           partial coherence analyses of fMRI data Sun 2004
function Cluster = estimate_cohe(Cluster,outdir,options,sel)
%load(fullfile(outdir,'parameter.mat'));
fprintf('Berechnung der Coherence\n');
idx = 0;
cmap = colormap(lines(size(Cluster,1)));
parameter = options.parameter;
coh_string = ['Coh' sel];
cohf_string = ['CohF' sel];
y_string = ['Y' sel '_mean'];
%%%%%%%%%%%%%%%%%%%%%%%
if options.Cluster2Voxel
    fprintf('coherence of Cluster2Voxel is not implemented yet\n');
end   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options.Cluster2Cluster
%    warning off
    if size(Cluster,1)>1
         for j=1:size(Cluster,1)
            x = Cluster{j,1}.(y_string);
            for k=1:size(Cluster,1)
                y = Cluster{k,1}.(y_string);
                %y = mean_tc_akt(:,k);
                [Cxy,F] = ck_rc3_get_coherence(x,...
                y,parameter.TR,parameter.hpf,parameter.lpf);
%                assignin('base','Cxy',Cxy);
%                assignin('base','F',F);

                Cluster{1,1}.ALL_C2C.(coh_string)(j,k,:)=Cxy;
                Cluster{1,1}.ALL_C2C.(cohf_string)(j,k,:)=F;
            end
            
            Cluster{j,1}.Cluster2Cluster.(coh_string) = squeeze(Cluster{1,1}.ALL_C2C.(coh_string)(j,:,:));
            Cluster{j,1}.Cluster2Cluster.(cohf_string) = squeeze(Cluster{1,1}.ALL_C2C.(cohf_string)(j,:,:));
         end        
    else
        msgbox('you must define at least 2 Clusters for Cluster2Cluster analysis');
    end
 %   warning on

    
    
    
end
if options.Voxel2Voxel
    fprintf('coherence of Voxel2Voxel is not implemented jet\n');
end
warning on
%assignin('base','R',R{idx});
%assignin('base','P',P{idx});
%assignin('base','Z',Z{idx});
%save(fullfile(outdir,['Cluster_' zname '.mat']),'Cluster','-v7.3');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% partial coherence analysis
function Cluster = estimate_part_cohe(Data,Cluster,V,outdir,options,zname)
load(fullfile(outdir,'parameter.mat'));
fprintf('Berechnung der partiellen Coherence\n');
idx = 0;
cmap = colormap(lines(size(Cluster,1)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options.Cluster2Voxel
    fprintf('coherence of Cluster2Voxel is not implemented jet\n');
end   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options.Cluster2Cluster
    warning off
    fprintf('schreibe Cluster');
    [Cluster_akt, Cluster_pas] = partial_Cluster(Cluster,outdir);
    assignin('base','Cluster_akt',Cluster_akt);
    assignin('base','Cluster_pas',Cluster_pas);
    if size(Cluster_akt,1)>1
        mean_tc = zeros(size(Cluster_akt{1,1}.Y_mean,1),size(Cluster_akt,1));
        for j=1:size(Cluster_akt,1)
            mean_tc(:,j) = Cluster_akt{j,1}.Y_mean;
        end
        for j=1:size(Cluster_akt,1)
            mean_tc_akt = mean_tc(:,[1:j-1 j+1:size(Cluster_akt,1)]);
            idx=idx+1;
            %%% nun die Coherenceanalyse 
            [Cxy,F] = ck_rc_get_coherence(Cluster_akt{j,1}.Y_mean,...
                mean_tc_akt,parameter.TR,0.01,0.2);
            %figure('Name',Cluster_akt{j,1}.name);
            %size(Cxy)
            %leg{1,1} = 'xxx';
            for kk=1:max(size(Cxy))
                dummy = smooth(Cxy{kk},20);
                plot(F{kk}(5:end-5),dummy(5:end-5),'Color',cmap(kk,:));
                hold on
                %leg = [leg; Cluster_akt{kk,1}.name];
            end
            leg = get_legend(j,Cluster_akt);
            legend(leg);
            hold off
        end
            % nun die passiven Cluster
            
            
                    mean_tc = zeros(size(Cluster_pas{1,1}.Y_mean,1),size(Cluster_pas,1));
        for j=1:size(Cluster_pas,1)
            mean_tc(:,j) = Cluster_pas{j,1}.Y_mean;
        end
        for j=1:size(Cluster_pas,1)
            mean_tc_akt = mean_tc(:,[1:j-1 j+1:size(Cluster_pas,1)]);
            idx=idx+1;
            %%% nun die Coherenceanalyse 
            [Cxy,F] = ck_rc_get_coherence(Cluster_pas{j,1}.Y_mean,...
                mean_tc_akt,parameter.TR,0.01,0.2);
           % figure('Name',Cluster_pas{j,1}.name);
            %size(Cxy)
            %leg{1,1} = 'xxx';
            for kk=1:max(size(Cxy))
                dummy = smooth(Cxy{kk},20);
                plot(F{kk}(5:end-5),dummy(5:end-5),'Color',cmap(kk,:));
                hold on
                %leg = [leg; Cluster_pas{kk,1}.name];
            end
            leg = get_legend(j,Cluster_pas);
            legend(leg);
            hold off
        end
    else
        msgbox('you must define at least 2 Clusters for Cluster2Cluster analysis');
    end
end
if options.Voxel2Voxel
    fprintf('coherence of Voxel2Voxel is not implemented jet\n');
end
warning on
%assignin('base','R',R{idx});
%assignin('base','P',P{idx});
%assignin('base','Z',Z{idx});
%save(fullfile(outdir,['Cluster_' zname '.mat']),'Cluster');
%save(fullfile(outdir,datafile),'Cluster','-v7.3');


function [T_string] = get_trans_string(y_string)

if strcmp(y_string,'Y_mean')
    T_string = 'T';
elseif strcmp(y_string,'Yc1_mean')
    T_string = 'Tc1';
elseif strcmp(y_string,'Yc2_mean')
    T_string = 'Tc2';
else
    fprintf('invalid string in estimate corr in ck_rc3_estimate_connectivity\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% TRANSINFORAMTION ANALYSIS %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Cluster = estimate_transinformation(Cluster,outdir,options,zname,sel)
%load(fullfile(outdir,'parameter.mat'));
fprintf('Berechnung der Transinformation\n');
idx = 0;
cmap = colormap(lines(size(Cluster,1)));
bins=options.Transinformation_bins;
lag=options.Transinformation_lag;
%parameter.bins=bins;
T_string = ['T' sel];
y_string = ['Y' sel '_mean'];

%%%%%%%%%%%%%%%%%%%%%%%
if options.Cluster2Voxel
    fprintf('transinformation of Cluster2Voxel is not implemented jet\n');
end   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options.Cluster2Cluster
    warning off
    if size(Cluster,1)>1
        mean_tc = zeros(length(Cluster{1,1}.(y_string)),length(Cluster));
        for j=1:size(Cluster,1)
            mean_tc(:,j) = Cluster{j,1}.(y_string);
        end
        for j=1:size(Cluster,1)
            mean_tc_akt = mean_tc(:,[1:j-1 j+1:size(Cluster,1)]);
            idx=idx+1;

            %%% nun die Transinformationsanalyse 
            %for k=1:size(mean_tc_akt,2)
            x = Cluster{j,1}.(y_string);
            for k=1:size(Cluster,1)
                y = Cluster{k,1}.(y_string);
                %y = mean_tc_akt(:,k);
                CMIF = ck_rc3_get_transinformation(x,y,lag,bins);

                %assignin('base','CMIF',CMIF);
                Cluster{1,1}.ALL_C2C.(T_string)(j,k,:)=CMIF(:,2);                
               
            end
            Cluster{j,1}.Cluster2Cluster.(T_string) = squeeze(Cluster{1,1}.ALL_C2C.(T_string)(j,:,:))';
            

        end
        
    else
        msgbox('you must define at least 2 Clusters for Cluster2Cluster analysis');
    end
    
end
if options.Voxel2Voxel
    fprintf('transinformation of Voxel2Voxel is not implemented jet\n');
end
warning on
%assignin('base','R',R{idx});
%assignin('base','P',P{idx});
%assignin('base','Cluster',Cluster);
%save(fullfile(outdir,['Cluster_transinformation_tmp.mat']),'Cluster');
%save(fullfile(outdir,['Cluster_' zname '.mat']),'Cluster','-v7.3');
%save(fullfile(outdir,'parameter.mat'),'parameter');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% SYNCHRONIZATION ANALYSIS %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%function Cluster = estimate_synchronization(Data,Cluster,V,outdir,options,zname)
function Cluster = estimate_synchronization(Cluster,outdir,options,zname,sel)
fprintf('Berechnung der Phasensynchronization mittels Kopplungsindex\n');
fprintf('Berechnung der Synchronization\n');

syn_string = ['Syn' sel];
plv_string = ['PLV' sel];
y_string = ['Y' sel '_mean'];
bins=options.Synchronization_bins;

%%%%%%%%%%%%%%%%%%%%%%%
if options.Cluster2Voxel
    fprintf('transinformation of Cluster2Voxel is not implemented jet\n');
end   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options.Cluster2Cluster
    %warning off
    if size(Cluster,1)>1
         for j=1:size(Cluster,1)
            x = Cluster{j,1}.(y_string);
            for k=1:size(Cluster,1)
                y = Cluster{k,1}.(y_string);
                %y = mean_tc_akt(:,k);
                [ci, plv] = ck_rc3_get_synchronization_ci_PLV(x,y,bins);
%                ci = ck_rc2_get_synchronization_kopplungsindex(x,y,bins);
                %assignin('base','CMIF',CMIF);
                Cluster{1,1}.ALL_C2C.(syn_string)(j,k)=ci;
                Cluster{1,1}.ALL_C2C.(plv_string)(j,k)=plv;
            end
            Cluster{j,1}.Cluster2Cluster.(syn_string) = squeeze(Cluster{1,1}.ALL_C2C.(syn_string)(j,:))';
            Cluster{j,1}.Cluster2Cluster.(plv_string) = squeeze(Cluster{1,1}.ALL_C2C.(plv_string)(j,:))';
         end        
    else
        msgbox('you must define at least 2 Clusters for Cluster2Cluster analysis');
    end
    
end
if options.Voxel2Voxel
    fprintf('transinformation of Voxel2Voxel is not implemented jet\n');
end
warning on
%assignin('base','R',R{idx});
%assignin('base','P',P{idx});
assignin('base','Cluster',Cluster);
%save(fullfile(outdir,['Cluster_transinformation_tmp.mat']),'Cluster');
%save(fullfile(outdir,['Cluster_' zname '.mat']),'Cluster','-v7.3');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% COMPLEXITY ANALYSIS %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%function Cluster = estimate_synchronization(Data,Cluster,V,outdir,options,zname)
function Cluster = estimate_complexity(Cluster,outdir,options,zname,sel)
%load(fullfile(outdir,'parameter.mat'));
fprintf('Berechnung der complexity mittels Kopplungsindex\n');

comp_string = ['Comp' sel];
y_string = ['Y' sel '_mean'];
%%%%%%%%%%%%%%%%%%%%%%%
if options.Cluster2Voxel
    fprintf('complextity of Cluster2Voxel is not implemented jet\n');
end   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options.Cluster2Cluster
    %warning off
    if size(Cluster,1)>1
         for j=1:size(Cluster,1)
            x = Cluster{j,1}.(y_string);
%             for k=1:size(Cluster,1)
%                 y = Cluster{k,1}.Y_mean;
                %y = mean_tc_akt(:,k);
                [entr] = ck_rc3_approx_entropy(options.Complexity_bins ,options.Complexity_r ,x);
%                ci = ck_rc2_get_synchronization_kopplungsindex(x,y,bins);
                %assignin('base','CMIF',CMIF);
%                Cluster{1,1}.ALL_C2C.Comp(j,k)=ci;
                
%            end
%            Cluster{j,1}.Cluster2Cluster.Syn = squeeze(Cluster{1,1}.ALL_C2C.Syn(j,:))';
            Cluster{j,1}.Cluster2Cluster.(comp_string) = entr;
         end        
    else
        msgbox('you must define at least 2 Clusters for Cluster2Cluster analysis');
    end
    
end
if options.Voxel2Voxel
    fprintf('complexity of Voxel2Voxel is not implemented jet\n');
end
%warning on
%assignin('base','R',R{idx});
%assignin('base','P',P{idx});
assignin('base','Cluster',Cluster);
%save(fullfile(outdir,['Cluster_transinformation_tmp.mat']),'Cluster');
%save(fullfile(outdir,['Cluster_' zname '.mat']),'Cluster','-v7.3');


function leg = get_legend(j,Cluster)
% Funktion extrahiert alle Clustername außer j
% und gibt diese als Spalten cell Array zurück
idx = 1;
for i=1:size(Cluster,1)
    if i~=j
        leg{idx,1}=Cluster{i,1}.name;
        idx = idx + 1;
    end
end

function [Cluster_akt, Cluster_pas] = partial_Cluster(Cluster,outdir)
%%%%%%%%%%%%%%%%%%%%%%%
%%% Aufspaltung der Daten nach aktiven und passiven Teilen 
%%%  mit anschließender Konkatenierung
load(fullfile(outdir,'parameter.mat'));
load(parameter.spmfile);
% Implementiert für einen Block der im SPM File als
% erster Kontrast eingegeben wurde
Cluster_akt = Cluster;
Cluster_pas = Cluster;
ons = SPM.Sess.U.ons;
dur = SPM.Sess.U.dur;
% erstelle einen Vektor mit 1 für aktiv und 0 für passiv
akt_vec = zeros(size(Cluster{1,1}.Y,1),1);
%for i=1:size(Cluster{1,1}.Y,1)
for j=1:size(ons,1)
    for k=1:dur(j,1)
        akt_vec(ons(j,1)+k,1)=1;
    end
end

for i=1:size(Cluster,1)
    Cluster_akt{i,1}.name = [Cluster{i,1}.name '_akt'];
    Cluster_pas{i,1}.name = [Cluster{i,1}.name '_pas'];
    Cluster_akt{i,1}.Y = zeros(1,size(Cluster{i,1}.Y,2));
    Cluster_pas{i,1}.Y = zeros(1,size(Cluster{i,1}.Y,2));
    Cluster_akt{i,1}.Y_mean = zeros(1,1);
    Cluster_pas{i,1}.Y_mean = zeros(1,1);
   
    idx_akt = 1;
    idx_pas = 1;
    for j=1:size(Cluster{i,1}.Y,1)
        if akt_vec(j,1)
            Cluster_akt{i,1}.Y(idx_akt,:) =  Cluster{i,1}.Y(j,:);
            Cluster_akt{i,1}.Y_mean(idx_akt,:) =  Cluster{i,1}.Y_mean(j,:);
            idx_akt = idx_akt + 1;
        else
            Cluster_pas{i,1}.Y(idx_pas,:) =  Cluster{i,1}.Y(j,:);
            Cluster_pas{i,1}.Y_mean(idx_pas,:) =  Cluster{i,1}.Y_mean(j,:);
            idx_pas = idx_pas + 1;
        end
    end
end


% return the identifier after the Cluster filename
function postfix = get_postfix(datafile)
[p,f,e]=fileparts(datafile);
t0 = findstr(f,'Cluster_');
if t0
    postfix = f(t0+8:end);
else
    error('Cluster files MUST be names as *Cluster_* and an identifier behind Cluster_');
end