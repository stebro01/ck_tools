function Network2 = ck_rc3_estimate_C2C_network_groups_statistic2(outdir,unc,...
    fwe,fdr, data, numgroup, method, normalize,C2C_stat,Cluster)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C2C netzwerk Gruppen statistic von einer oder von 2 Gruppen
%  Berechnet werden Unterschiede in vorher definierten Netzwerken
% übergeben werden
% outdir ...        das aktuelle Verzeichnis in dem gearbeitet wird
% unc,fdr,fwe ...   die Korrekturschwellen für die Statistik
% data...           Datenstruktur mit den data.R1, data.Z1... 
%                                         data.R2, data.Z2
%                   Z hat 3 dim (numregionen, numregionen, numsubj)
% numgroup...       die Anzahl der Gruppen ( 1 oder 2)
% C2C_stat...       die bisher bereits errechnete Gruppenstatistik
% output
% C2C_stat ...      wobei im Wesentlichen die Variable Network
%                   bearbeitet wurde


if normalize 
   data = normalize_data(data);
end


 
str = '';
if method == 3
    str = sprintf('Network statistic paired ttest %i vs %i; p = %g',size(data.Z1,3),size(data.Z2,3),unc);
elseif method == 2
    str = sprintf('Network statistic 2s ttest %i vs %i; p = %g',size(data.Z1,3),size(data.Z2,3),unc);
elseif method == 1
    str = sprintf('Network statistic 1s ttest %i ; p = %g',size(data.Z1,3),unc);
end

 h=waitbar(0,['Please wait...' str]);

 
% generiere Gruppen Verbindungsmatrix der Netzwerke GNM
GNM = get_network_matrix(Cluster);
% generiere Gruppen Verbindungsmatrix zwischen den Netzwerken
[GiNM,inter_network_conn] = get_inter_network_matrix(Cluster);

%assignin('base','RGiN1',RGiN1);
%assignin('base','RGiN2',RGiN2);

Network2.name = sprintf('2s intra Netzwerk ttest %i vs %i; p = %g',size(data.Z1,3),size(data.Z2,3),unc);
Network2.networknames = Cluster{1,1}.net.global.network_names;
[Network2.network_indices] = get_network_names_numbers(Cluster);
num_nets = length(Network2.networknames);
Network2.RGN1 = data.RN1;
Network2.RGN2 = data.RN2;
Network2.ZGN1 = data.ZN1;
Network2.ZGN2 = data.ZN2;
Network2.RGiN1 = data.RiN1;
Network2.RGiN2 = data.RiN2;
Network2.ZGiN1 = data.ZiN1;
Network2.ZGiN2 = data.ZiN2;


if numgroup==2
    % Statistik fuer Netzwerke
    
    % two- sample oder paired t-statistic
    Network2= get_Z_statistic(Network2,unc,method);
    
    
end
close(h);




function Network2 = get_Z_statistic(Network2,unc,method)


for i=1:size(Network2.ZGN1,1) % Networks
    Network2.MeanZdiff(i)=mean(Network2.ZGN1(i,:))-mean(Network2.ZGN2(i,:));
    Network2.MeanRdiff(i)=mean(Network2.RGN1(i,:))-mean(Network2.RGN2(i,:));
    % einfache Netzwerke / Intranetwork
    if method == 2    % two sample t-test
        
        [Network2.Zh(i), Network2.Zp(i), Network2.Zci(i,:),...
            Zstats] = ttest2(Network2.ZGN1(i,:),Network2.ZGN2(i,:),unc,'both');
        Network2.Zt(i)=Zstats.tstat;
        Network2.Zsd(i)=Zstats.sd;
        
        [Network2.Z1k2h(i), Network2.Z1k2p(i)]=ttest2(Network2.ZGN1(i,:),...
            Network2.ZGN2(i,:),unc,'left');
        [Network2.Z2k1h(i), Network2.Z2k1p(i)]=ttest2(Network2.ZGN1(i,:),...
            Network2.ZGN2(i,:),unc,'right');
        
    end
    if method == 3   % paired t-test
        [Network2.Zh(i), Network2.Zp(i), Network2.Zci(i,:),...
            Zstats] = ttest(Network2.ZGN1(i,:),Network2.ZGN2(i,:),unc,'both');
        Network2.Zt(i)=Zstats.tstat;
        Network2.Zsd(i)=Zstats.sd;
        [Network2.Z1k2h(i), Network2.Z1k2p(i)]=ttest(Network2.ZGN1(i,:),...
            Network2.ZGN2(i,:),unc,'left');
        [Network2.Z2k1h(i), Network2.Z2k1p(i)]=ttest(Network2.ZGN1(i,:),...
            Network2.ZGN2(i,:),unc,'right');
        
    end
    % nun noch internetwork
    for j = 1:size(Network2.ZGiN1,2)
        if method ==2

            [Network2.ZiNh(i,j), Network2.ZiNp(i,j), Network2.ZiNci(i,j,:),...
                ZiNstats] = ttest2(squeeze(Network2.ZGiN1(i,j,:)),squeeze(Network2.ZGiN2(i,j,:)),unc,'both');
        end
        if method == 3
            [Network2.ZiNh(i,j), Network2.ZiNp(i,j), Network2.ZiNci(i,j,:),...
                ZiNstats] = ttest(Network2.ZGiN1(i,j,:),Network2.ZGiN2(i,j,:),unc,'both');
        end
        Network2.ZiNt(i,j)=Zstats.tstat;
        Network2.ZiNsd(i,j)=Zstats.sd;
        
        Network2.MeanZiNdiffq(i,j)=mean(Network2.ZGiN1(i,j,:))-mean(Network2.ZGiN2(i,j,:));
        Network2.MeanRiNdiffq(i,j)=mean(Network2.RGiN1(i,j,:))-mean(Network2.RGiN2(i,j,:));
    end
    for kk=1:size(Network2.MeanZiNdiffq,1)
            Network2.MeanZiNdiffq(kk,kk) = Network2.MeanZdiff(kk);
            Network2.MeanRiNdiffq(kk,kk) = Network2.MeanRdiff(kk);
            Network2.ZiNp(kk,kk)=Network2.Zp(kk);
    end
end





function [Rq] = get_internetwork_square(R,conn,num_nets)
Rq = zeros(num_nets,num_nets,size(R,1));
for j = 1:size(conn,1) % internetwork konnnections
    for i=1:size(R,1) % Subjects
        Rq(conn(j,1),conn(j,2),i) = R(i,j);
        Rq(conn(j,2),conn(j,1),i) = R(i,j);
    end
end

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


function G = get_network_matrix(Cluster)

num_global_networks = Cluster{1}.net.global.num_global_networks;
N_conn = Cluster{1}.net.global.N_conn;
N_conn = N_conn(:,1:num_global_networks);
%    g(i)=Cluster{i}.group;

num_c = length(Cluster);
G = zeros(num_c,num_c,num_global_networks);
%num_global_networks
for k=1:num_global_networks
for i=1:num_c
    for j=1:num_c

        %assignin('base','N_conn',N_conn)
        %assignin('base','G',G)
        
        if i~=j && N_conn(i,k) && N_conn(j,k)
            G(i,j,k)=1;
        end
    end
end
end

% Matrix von Verbindenung zwischen Netzwerken
function [G,inter_network_conn] = get_inter_network_matrix(Cluster)
% G ist die verbindungsmatrix mit Dimension cluster x Cluster x
% inter_network_index
% inter_network_conn ... welches Netzwerk mit welchem in der Matrix
% G(:,:,k) steht 
% Bsp.: G(1,3,14)=1 inter_network_conn(14,:)=[5,7] .... Der Cluster 1 aus
% Netzwerk 5 hat eine Verbindung mit Cluster 3 aus Netzwerk 7, eingetragen
% in der 14 Matrix von G
% 

num_global_networks = Cluster{1}.net.global.num_global_networks;
N_conn = Cluster{1}.net.global.N_conn;
N_conn = N_conn(:,1:num_global_networks);
%    g(i)=Cluster{i}.group;

num_c = length(Cluster);
G = zeros(num_c,num_c,num_global_networks);
G = zeros(size(Cluster,1));
inter_network_index = 0;
for k=1:num_global_networks
    for l=k+1:num_global_networks
        inter_network_index = inter_network_index + 1;
        inter_network_conn(inter_network_index,1:2)=[k,l];
        for i=1:num_c
            for j=1:num_c
            
            if i~=j && N_conn(i,k) && N_conn(j,l)
                G(i,j,inter_network_index)=1;
            end
        end
    end
end
end




% mittle die R-Werte fuer jedes Netzwerk fuer jeden probanden
function [RG RGstd] = get_group_mean_R(R,GM)

RG = zeros(size(R,3),size(GM,3));

for p=1:size(R,3) % ueber die Probanden
    for g=1:size(GM,3) % ueber die Netzwerke
        dummysum = 0;
        dummyidx = 0;
        X=0;
        for i=1:size(GM,1) % ueber die Netzwerkmatrix
            for j=1:size(GM,2) % ueber die Netzwerkmatrix
                if GM(i,j,g)==1
                    %fprintf('%.0f  %.0f  %.0f \n',i,j,p)
                    dummysum = dummysum + R(i,j,p);
                    dummyidx = dummyidx + 1;
                    X(dummyidx)=R(i,j,p);
                end
            end
        end
        RG(p,g)= mean(X);
        RGstd(p,g) = std(X);
%         str = '';
%         str = sprintf('group %.0f',g);
%         names{g+1} = str;
    end
end




function [networkindices] = get_network_names_numbers(Cluster)

networknumbers = [1:length(Cluster{1,1}.net.global.num_global_networks)];
N = Cluster{1,1}.net.global.N_conn;
for i=1:Cluster{1,1}.net.global.num_global_networks
idx = 1;
for j=1:length(Cluster)
    if N(j,i)
        networkindices{i}(idx)= j;
        idx = idx + 1;
    end
end
end
    
% % network indices
% idx = 1;
% indi = ones(1,length(networknumbers));
% for i=1:size(Cluster,1)
%     x=Cluster{i,1}.net.local.group; % Cluster 1 gehoert zu Netzwerk 5
%     networkindices{x}(indi(x))=i; % networkindices{5}(5) = clusternumber
%     indi(x) = indi(x)+1;
% end


% 
% function [networknames, networknumbers, networkindices] = get_network_names_numbers(Cluster)
% 
% for i=1:size(Cluster,1)
%     g(i)=Cluster{i,1}.net.local.group;
%     N{i}=Cluster{i,1}.net.loccal.groupname;
% 
% end
% 
% idx = 1;
% for i=1:max(g)
%     suche = 1;
%     for j=1:length(g)
%         if g(j)==i && suche
%             networknumbers(idx)=g(j);
%             networknames{idx} = N{j};
%             idx = idx +1;
%             suche = 0;
%         end
%     end
% end
% 
% % network indices
% idx = 1;
% indi = ones(1,length(networknumbers));
% for i=1:size(Cluster,1)
%     x=Cluster{i,1}.net.local.group;
%     networkindices{x}(indi(x))=i;
%     indi(x) = indi(x)+1;
% end
% 


