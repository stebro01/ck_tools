function C2C_stat = ck_rc3_estimate_C2C_statistic(outdir,unc,...
    fwe,fdr, data, numgroup, method, normalize)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C2C statistic von einer oder von 2 Gruppen
% übergeben werden
% outdir ...        das aktuelle Verzeichnis in dem gearbeitet wird
% unc,fdr,fwe ...   die Korrekturschwellen für die Statistik
% data...           Datenstruktur mit den data.R1, data.Z1... 
%                                         data.R2, data.Z2
%                   Z hat 3 dim (numregionen, numregionen, numsubj)
% numgroup...       die Anzahl der Gruppen ( 1 oder 2)

%%%%%%%%%%%%%% aus early blind paper
%The individual z scores were entered into a one-sample two-tailed t-test 
%to determine if the two brain regions show significant functional 
%connectivity within each group. They were also entered into a 
%two-sample two-tailed t-test to determine if the functional 
%connectivities were significantly different between the two groups.
%A t-test was performed for all the 6670 (116 x 115/2) functional 
% connectivities, so a correction for multiple comparisons was strictly 
% necessary. The false discovery rate (FDR) approach was applied to find 
% a threshold that would restrict the expected proportion of type I errors 
% to lower than 0.05 (Benjamini and Yekutieli, 2001Go; Salvador et al., 
% 2005aGo). In this study, we identified the significant differences in 
% functional connectivities between the blind and sighted subjects 
% according to the following two criteria: (a) the z values were 
% significantly different from zero at least in one group at P < 0.05 
% (one-sample two-tailed t-test; FDR corrected); (b) the z scores were 
% significantly different between the two groups at P < 0.05 (two-sample 
% two-tailed t-test; FDR corrected).

if normalize 
   data = normalize_data(data);
end

dim1= size(data.Z1,1);
dim2= size(data.Z1,2);
C2C_stat.Zh=zeros(dim1,dim2);
C2C_stat.Zp=zeros(dim1,dim2);
C2C_stat.Zci=zeros(dim1,dim2,2);
C2C_stat.T=zeros(dim1,dim2);
C2C_stat.S=zeros(dim1,dim2);
C2C_stat.R1mean=zeros(dim1,dim2);
C2C_stat.R1std=zeros(dim1,dim2);
C2C_stat.R2mean=zeros(dim1,dim2);
C2C_stat.R2std=zeros(dim1,dim2);
C2C_stat.Syn1mean=zeros(dim1,dim2);
C2C_stat.Syn1std=zeros(dim1,dim2);
C2C_stat.Syn2mean=zeros(dim1,dim2);
C2C_stat.Syn2std=zeros(dim1,dim2);
C2C_stat.PLV1mean=zeros(dim1,dim2);
C2C_stat.PLV1std=zeros(dim1,dim2);
C2C_stat.PLV2mean=zeros(dim1,dim2);
C2C_stat.PLV2std=zeros(dim1,dim2);
 
str = '';
if method == 3
    str = sprintf('paired ttest %i vs %i; p = %g',size(data.Z1,3),size(data.Z2,3),unc);
elseif method == 2
    str = sprintf('2s ttest %i vs %i; p = %g',size(data.Z1,3),size(data.Z2,3),unc);
elseif method == 1
    str = sprintf('1s ttest %i ; p = %g',size(data.Z1,3),unc);
end

 h=waitbar(0,['Please wait...' str]);  

if numgroup==2
    for i=1:dim1
        waitbar(i/dim1)
        for j=1:size(data.Z1,2)
%             assignin('base','data',data);
%             squeeze(data.Z1(i,j,:))
%             squeeze(data.Z2(i,j,:))
%             ttest2(squeeze(data.Z1(i,j,:)),...
%                 squeeze(data.Z2(i,j,:)),unc,'both')
%             
            if method == 2          %%% added by sb 280712
                C2C_stat.name = sprintf('2s ttest %i vs %i; p = %g',size(data.Z1,3),size(data.Z2,3),unc);
                [C2C_stat.Zh(i,j), C2C_stat.Zp(i,j), C2C_stat.Zci(i,j,:),...
                    Zstats] = ttest2(squeeze(data.Z1(i,j,:)),...
                    squeeze(data.Z2(i,j,:)),unc,'both');
                C2C_stat.Zt(i,j)=Zstats.tstat;
                C2C_stat.Zsd(i,j)=Zstats.sd;

                [C2C_stat.Z1k2h(i,j), C2C_stat.Z1k2p(i,j)]=ttest2(squeeze(data.Z1(i,j,:)),...
                    squeeze(data.Z2(i,j,:)),unc,'left');
                [C2C_stat.Z2k1h(i,j), C2C_stat.Z2k1p(i,j)]=ttest2(squeeze(data.Z1(i,j,:)),...
                    squeeze(data.Z2(i,j,:)),unc,'right');
                C2C_stat.MeanZdiff(i,j)=mean(squeeze(data.Z1(i,j,:)))-mean(squeeze(data.Z2(i,j,:)));
                C2C_stat.stdg1(i,j)=std(squeeze(data.Z1(i,j,:)));
                C2C_stat.stdg2(i,j)=std(squeeze(data.Z2(i,j,:)));

                C2C_stat.R1mean(i,j)=mean(squeeze(data.R1(i,j,:)));
                C2C_stat.R1std(i,j)=std(squeeze(data.R1(i,j,:)));
                C2C_stat.R2mean(i,j)=mean(squeeze(data.R2(i,j,:)));
                C2C_stat.R2std(i,j)=std(squeeze(data.R2(i,j,:)));
                C2C_stat.Syn1mean(i,j)=mean(squeeze(data.Syn1(i,j,:)));
                C2C_stat.Syn1std(i,j)=std(squeeze(data.Syn1(i,j,:)));
                C2C_stat.Syn2mean(i,j)=mean(squeeze(data.Syn2(i,j,:)));
                C2C_stat.Syn2std(i,j)=std(squeeze(data.Syn2(i,j,:)));
                C2C_stat.MeanSyndiff(i,j)=mean(squeeze(data.Syn1(i,j,:)))-mean(squeeze(data.Syn2(i,j,:)));
                [C2C_stat.Synh(i,j), C2C_stat.Synp(i,j), C2C_stat.Synci(i,j,:),...
                    C2C_stat.Synstats] = ttest2(squeeze(data.Syn1(i,j,:)),...
                    squeeze(data.Syn2(i,j,:)),unc,'both');

                C2C_stat.PLV1mean(i,j)=mean(squeeze(data.PLV1(i,j,:)));
                C2C_stat.PLV1std(i,j)=std(squeeze(data.PLV1(i,j,:)));
                C2C_stat.PLV2mean(i,j)=mean(squeeze(data.PLV2(i,j,:)));
                C2C_stat.PLV2std(i,j)=std(squeeze(data.PLV2(i,j,:)));
                C2C_stat.MeanPLVdiff(i,j)=mean(squeeze(data.PLV1(i,j,:)))-mean(squeeze(data.PLV2(i,j,:)));
                [C2C_stat.PLVh(i,j), C2C_stat.PLVp(i,j), C2C_stat.PLVci(i,j,:),...
                    C2C_stat.PLVstats] = ttest2(squeeze(data.PLV1(i,j,:)),...
                    squeeze(data.PLV2(i,j,:)),unc,'both');

                
                C2C_stat.data = data;

            end%%%               added by sb 280712
            
            %%% added by sb 280712
            if method==3
                % paired test
                C2C_stat.name = sprintf('paired ttest %i vs %i; p = %g',size(data.Z1,3),size(data.Z2,3),unc);

                [C2C_stat.Zh(i,j), C2C_stat.Zp(i,j), C2C_stat.Zci(i,j,:),...
                    Zstats] = ttest(squeeze(data.Z1(i,j,:)),...
                    squeeze(data.Z2(i,j,:)),unc,'both');
                C2C_stat.Zt(i,j)=Zstats.tstat;
                C2C_stat.Zsd(i,j)=Zstats.sd;

                [C2C_stat.Z1k2h(i,j), C2C_stat.Z1k2p(i,j)]=ttest(squeeze(data.Z1(i,j,:)),...
                    squeeze(data.Z2(i,j,:)),unc,'left');
                [C2C_stat.Z2k1h(i,j), C2C_stat.Z2k1p(i,j)]=ttest(squeeze(data.Z1(i,j,:)),...
                    squeeze(data.Z2(i,j,:)),unc,'right');
                C2C_stat.MeanZdiff(i,j)=mean(squeeze(data.Z1(i,j,:)))-mean(squeeze(data.Z2(i,j,:)));
                C2C_stat.stdg1(i,j)=std(squeeze(data.Z1(i,j,:)));
                C2C_stat.stdg2(i,j)=std(squeeze(data.Z2(i,j,:)));

                C2C_stat.R1mean(i,j)=mean(squeeze(data.R1(i,j,:)));
                C2C_stat.R1std(i,j)=std(squeeze(data.R1(i,j,:)));
                C2C_stat.R2mean(i,j)=mean(squeeze(data.R2(i,j,:)));
                C2C_stat.R2std(i,j)=std(squeeze(data.R2(i,j,:)));
                C2C_stat.data = data;
            end
            %%% END added by sb 280712
            
        end
    end
end


if numgroup==1
    for i=1:dim1
        waitbar(i/dim1)
        for j=1:size(data.Z1,2)
              assignin('base','data',data);
              
%             squeeze(data.Z1(i,j,:))
%             squeeze(data.Z2(i,j,:))
%             ttest2(squeeze(data.Z1(i,j,:)),...
%                 squeeze(data.Z2(i,j,:)),unc,'both')
%             
            [C2C_stat.Zh(i,j), C2C_stat.Zp(i,j), C2C_stat.Zci(i,j,:),...
                Zstats] = ttest(squeeze(data.Z1(i,j,:)),0,unc,'both');
            C2C_stat.Zt(i,j)=Zstats.tstat;
            C2C_stat.Zsd(i,j)=Zstats.sd;
            

            
            
            C2C_stat.R_std(i,j)=std(squeeze(data.R1(i,j,:)));
            C2C_stat.T_std(i,j)=std(squeeze(data.T1(i,j,:)));
            
            C2C_stat.S_std(i,j)=std(squeeze(data.S1(i,j,:)));
            C2C_stat.data = data;
            
            
%             [C2C_stat.Z1k2h(i,j), C2C_stat.Z1k2p(i,j)]=ttest2(squeeze(data.Z1(i,j,:)),...
%                 squeeze(data.Z2(i,j,:)),unc,'left');
%             [C2C_stat.Z2k1h(i,j), C2C_stat.Z2k1p(i,j)]=ttest2(squeeze(data.Z1(i,j,:)),...
%                 squeeze(data.Z2(i,j,:)),unc,'right');
%             C2C_stat.MeanZdiff(i,j)=mean(squeeze(data.Z1(i,j,:)))-mean(squeeze(data.Z2(i,j,:)));
%             C2C_stat.stdg1(i,j)=std(squeeze(data.Z1(i,j,:)));
%             C2C_stat.stdg2(i,j)=std(squeeze(data.Z2(i,j,:)));

        end
        
    end
                C2C_stat.R_mean=mean(data.R1,3);
            C2C_stat.T_mean=mean(data.T1,3);
            C2C_stat.S_mean=mean(data.S1,3);
            C2C_stat.S_mean=mean(data.S1,3);
end



close(h)


function data = normalize_data(data)
% Die Datensätze werden um den Absolutbetrag der summe aller R-Werte 
% korrigiert um z.B. altersabhängige Steigerungen oder Senkungen der 
% Connektivität die lokal unspezifisch sind zu kompensieren
%  dann könnte ich aber keine R-Z-Transformation mehr durchführen
% also korrigiere ich einfach nur die Z-Werte
%data=data;
% ermittle den durchschnittlichen Wert der Z-Werte

    Z1=data.Z1;
    Z2=data.Z2;

%      assignin('base','Z1',Z1)
%      assignin('base','Z2',Z2)

%%%%%SB ADDED %%%%%%%%%%%%
%remove infinities >> aus irgendeinem Grund stehen in den Diagonalen INF
%und sehr große Werte >> die werden genullt, dann läuft die
%Normalize-Routine
    for i=1:size(Z1,3)
        for j = 1:size(Z1,2)
             Z1(j,j,i) = 0;
             Z2(j,j,i) = 0;
             Z1(Z1 == Inf) = 0;
             Z2(Z2 == Inf) = 0;
       end
    end 
%%%%END
    
    for i=1:size(Z1,3)
        x(i)=sum(sum(abs(Z1(:,:,i)))); 
    end
    mr=mean(x)/(size(Z1,1)*(size(Z1,2)-1));
    m=mean(x);

    for i=1:size(Z1,3)
        X1=squeeze(Z1(:,:,i));
        mi=sum(sum(abs(X1)));
        X1=X1*(m/mi);
        Z1(:,:,i)=X1;
    end
    for i=1:size(Z2,3)
        X2=squeeze(Z2(:,:,i));
        mi=sum(sum(abs(X2)));
        X2=X2*(m/mi);
        Z2(:,:,i)=X2;
    end
%      assignin('base','Z1b',Z1)
%      assignin('base','Z2b',Z2)
    data.Z1=Z1;
    data.Z2=Z2;



% Es fehlt noch eine Auswertung für die Gruppenstatistik von 
% der Transinformation und der Synchronization
% solange werden die Werte vorerst unverändert übernommne
% hieraus kann jedoch keine Statistik erstellt werden
% Eine Methode für die Transinformation wird in folgendem 
% paper präsentiert: Causal visual interactions as revealed by an information
% theoretic measure and fMRI


% 
% 
% if nargin<1
%     outdir = uigetdir();
% end
% if nargin<2
%     p_unc = 0.001;
% end
% if nargin<3
%     threshold_fwe = 0.05;
% end
% if nargin<4
%     threshold_fdr = 0.05;
% end
% 

% 
% c = 0;
% while true
% 
%     [f,p] = uigetfile('Cluster*.mat','Select Cluster-Files','MultiSelect','on');
%     if isnumeric(f)
%         break
%     end
%     c = c + 1;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % init der Erbebnistabelle
%     if size(f,2)>0
%         load(fullfile(p,f{1,1}));
%     else
%         load(fullfile(p,f));
%     end
%     C2C_stat{c,1} = struct; %#ok<AGROW>
%     C2C_stat{c,1}.R = zeros(size(Cluster,1),size(Cluster,1),size(f,2)); 
%     C2C_stat{c,1}.P = zeros(size(Cluster,1),size(Cluster,1),size(f,2)); %#ok<AGROW>
%     C2C_stat{c,1}.Z = zeros(size(Cluster,1),size(Cluster,1),size(f,2)); %#ok<AGROW>
%     for i=1:size(Cluster,1)
%         C2C_stat{c,1}.clusternames{i,1}=Cluster{i,1}.name; %#ok<AGROW>
%         C2C_stat{c,1}.clusterXYZ{i,1}=Cluster{i,1}.XYZ; %#ok<AGROW>
%     end
%     tmp = inputdlg('Cluster description');
%     C2C_stat{c,1}.name = tmp{1,1};
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%     % Schleife über alle ausgewählten Clusterstrukturen
%     % Jeder Proband hat eine Clusterstruktur
%     for i=1:size(f,2)
%         load(fullfile(p,f{1,i}));
%         %Schleife über die die einzelnen Cluster des Probanden
%         for j=1:size(Cluster,1)
%             C2C_stat{c,1}.R([1:j-1, j+1:size(Cluster,1)],j,i)= ...
%                 sym(Cluster{j,1}.Cluster2Cluster.R); %#ok<AGROW>
%             C2C_stat{c,1}.R(j,j,i)=1; %#ok<AGROW>
%             C2C_stat{c,1}.P([1:j-1, j+1:size(Cluster,1)],j,i)= ...
%                 Cluster{j,1}.Cluster2Cluster.P; %#ok<AGROW>
%             C2C_stat{c,1}.P(j,j,i)=1; %#ok<AGROW>
%             C2C_stat{c,1}.Z([1:j-1, j+1:size(Cluster,1)],j,i)= ...
%                 Cluster{j,1}.Cluster2Cluster.Z; %#ok<AGROW>
%             C2C_stat{c,1}.Z(j,j,i)=1; %#ok<AGROW>
%         end
%     end
% end
% 
% %%% Berechne die Statistic
% % one sample (two tailed) t-test für jede C2C_stat cell
% % und two sample t-test untereinander
% if nargin<5
%     n = (size(Cluster,1)-1) * (size(Cluster,1)-1);
% end
% p_fwe = threshold_fwe / n;
% for c=1:size(C2C_stat,1)
%     % one sample t-test
%     for i=1:size(C2C_stat{c,1}.Z,1)
%         for j=1:size(C2C_stat{c,1}.Z,2)
%             [h,p,ci,stats] = ttest(squeeze(C2C_stat{c,1}.Z(i,j,:)));
%             %C2C_stat{c,1}.one_sample_t_test.issignificant(i,j)=h;
%             C2C_stat{c,1}.one_sample_t_test.p(i,j)=p;
%             C2C_stat{c,1}.one_sample_t_test.ci1(i,j)=ci(1,1);
%             C2C_stat{c,1}.one_sample_t_test.ci2(i,j)=ci(2,1);
%             C2C_stat{c,1}.one_sample_t_test.stats{i,j}=stats;
%             C2C_stat{c,1}.one_sample_t_test.p_bonf(i,j)=p*(size(Cluster,1)-1);
%            % C2C_stat{c,1}.one_sample_t_test.issignificant_bonf(i,j)=...
%             %    p<(threshold/n);
%                 %(p*(size(Cluster,1)-1))<0.05;
%             % two sample t-test
%             if size(C2C_stat,1)>1
%                 for k=1:size(C2C_stat,1)
%                     [h,p,ci,stats] = ttest2(squeeze(C2C_stat{c,1}.Z(i,j,:)),...
%                         squeeze(C2C_stat{k,1}.Z(i,j,:)));
%                     %C2C_stat{c,1}.two_sample_t_test{k,1}.issignificant(i,j)=h;
%                     C2C_stat{c,1}.two_sample_t_test{k,1}.p(i,j)=p;
%                     C2C_stat{c,1}.two_sample_t_test{k,1}.ci1(i,j)=ci(1,1);
%                     C2C_stat{c,1}.two_sample_t_test{k,1}.ci2(i,j)=ci(2,1);
%                     C2C_stat{c,1}.two_sample_t_test{k,1}.stats{i,j}=stats;
%                     C2C_stat{c,1}.two_sample_t_test{k,1}.name= ...
%                         [C2C_stat{c,1}.name ' vs ' C2C_stat{k,1}.name];
%                     %C2C_stat{c,1}.two_sample_t_test{k,1}.p_bonf(i,j)= ...
%                      %   p*(size(Cluster,1)-1);
%                     %C2C_stat{c,1}.two_sample_t_test{k,1}.issignificant_bonf(i,j)= ...
%                      % p<(threshold/n);   
%                 %    p<((1-((1-threshold)^n))/n);
%                 end
%             end
%         end
%     end
%     % FWE
%     C2C_stat{c,1}.one_sample_t_test.issignificant_fwe=...
%         C2C_stat{c,1}.one_sample_t_test.p<p_fwe;
%     % FDR
%     p_fdr = get_fdr_threshold(C2C_stat{c,1}.one_sample_t_test.p,...
%                 n,threshold_fdr);
%     C2C_stat{c,1}.one_sample_t_test.issignificant_fdr=...
%         C2C_stat{c,1}.one_sample_t_test.p<p_fdr;
%     % uncorrected
%     C2C_stat{c,1}.one_sample_t_test.issignificant_unc=...
%         C2C_stat{c,1}.one_sample_t_test.p<p_unc;
% 
%     if size(C2C_stat,1)>1
%         for k=1:size(C2C_stat,1)
%             % FWE
%             C2C_stat{c,1}.two_sample_t_test{k,1}.issignificant_fwe=...
%                 C2C_stat{c,1}.two_sample_t_test{k,1}.p<p_fwe;
%             % FDR
%             p_fdr = get_fdr_threshold(C2C_stat{c,1}.two_sample_t_test{k,1}.p,...
%                 n,threshold_fdr);
%             C2C_stat{c,1}.two_sample_t_test{k,1}.issignificant_fdr=...
%                 C2C_stat{c,1}.two_sample_t_test{k,1}.p<p_fdr;
%             % uncorrected
%             C2C_stat{c,1}.two_sample_t_test{k,1}.issignificant_unc=...
%                 C2C_stat{c,1}.two_sample_t_test{k,1}.p<p_unc;
%         end
%     end
% 
% end
% 
% assignin('base','C2C_stat',C2C_stat);
% save(fullfile(outdir,'C2C_stat.mat'),'C2C_stat');


function p_fdr = get_fdr_threshold(P,n,threshold_fdr)
% berechne die False Discovery Rate
% sortiere die p Werte beginnend mit dem geringsten
% suche größtes p für das gilt: p < idx/n*threshold
% alle kleineren p-Werte sind signifikant
% 1. suche die Grenze
idx = 0;
for ii=1:size(P,1)
    for jj=1:size(P,2)
        if ii~=jj
            idx = idx + 1;
            p_list(idx) =  P(ii,jj);
        end
    end
end
p_list = sort(p_list);
gefunden = 0;
p_fdr = 0;
for ii=1:size(p_list,2)
    if p_list(ii)<=(ii/n*threshold_fdr)
        gefunden =1;
    else
        if gefunden
            gefunden = 0;
            p_fdr = p_list(ii-1);
        end
    end
end
%assignin('base','p_list',p_list);


