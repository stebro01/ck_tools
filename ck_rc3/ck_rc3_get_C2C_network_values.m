function C2C_stat_g1 = ck_rc3_get_C2C_network_values(outdir,C2C_stat_g1,Cluster)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C2C netzwerk Gruppen statistic von einer oder von 2 Gruppen
%  Berechnet werden Unterschiede in vorher definierten Netzwerken
% übergeben werden
% outdir ...        das aktuelle Verzeichnis in dem gearbeitet wird
% C2C_stat_g1...       die Values (R,Z, Syn,Power ... ) der Gruppe x


% generiere Gruppen Verbindungsmatrix der Netzwerke GNM
C2C_stat_g1.GNM = get_network_matrix(Cluster);
% generiere Gruppen Verbindungsmatrix zwischen den Netzwerken
[C2C_stat_g1.GiNM,C2C_stat_g1.inter_network_conn] = get_inter_network_matrix(Cluster);
C2C_stat_g1.networknames = Cluster{1,1}.net.global.network_names;
C2C_stat_g1.networknumbers = [1:length(C2C_stat_g1.networknames)];
[C2C_stat_g1.network_indices] = get_network_names_numbers(Cluster);

% mittle die R-Werte ueber jedes Netzwerk fuer jeden
% Probanden
%assignin('base','R1',data.R1);
[C2C_stat_g1.RN] = get_network_data(C2C_stat_g1.R,C2C_stat_g1.GNM,length(C2C_stat_g1.networknames));

[C2C_stat_g1.RiN, C2C_stat_g1.RiNstd] = get_internetwork_data(C2C_stat_g1.R,C2C_stat_g1.GiNM,length(C2C_stat_g1.networknames),C2C_stat_g1.inter_network_conn);


C2C_stat_g1.ZN = ck_rc3_fisher_transformation(C2C_stat_g1.RN);
C2C_stat_g1.ZiN = ck_rc3_fisher_transformation(C2C_stat_g1.RiN);
[C2C_stat_g1.GiNMq] = get_internetwork_square_GiNM(C2C_stat_g1.GiNM,C2C_stat_g1.inter_network_conn,length(C2C_stat_g1.networknames));



function [Mq] = get_internetwork_square_GiNM(M,conn,num_nets)
Mq = zeros(size(M,1),size(M,2),num_nets,num_nets);
for j = 1:size(conn,1) % internetwork konnnections
        Mq(:,:,conn(j,1),conn(j,2)) = M(:,:,j);
        Mq(:,:,conn(j,2),conn(j,1)) = Mq(:,:,j);
end


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

function [RG, RGstd] = get_network_data(R,GM,num_nets)

RG = zeros(num_nets,size(R,3));

RGstd = zeros(num_nets,size(R,3));

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
        RG(g,p)= mean(X);
        RGstd(g,p) = std(X);
    end
end




function [RG, RGstd] = get_internetwork_data(R,GiNM,num_nets,conn)

RG = zeros(num_nets,num_nets,size(R,3));
RGstd = zeros(num_nets,num_nets,size(R,3));
% assignin('base','R',R);
% assignin('base','GiNM',GiNM);
% assignin('base','num_nets',num_nets);
% assignin('base','conn',conn);
for p=1:size(R,3) % ueber die Probanden
    for iN = 1:size(GiNM,3)
        x = conn(iN,1);
        y = conn(iN,2);
        idx = squeeze(GiNM(:,:,iN));
        Rtmp = squeeze(R(:,:,p));
        RG(x,y,p) = mean(Rtmp(idx==1));
        RGstd(x,y,p) = std(Rtmp(idx==1));
        RG(y,x,p) = mean(Rtmp(idx==1));
        RGstd(y,x,p) = std(Rtmp(idx==1));
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
    


