function NET = ck_rc3_get_C2C_network_values3(outdir,C2C_stat,Cluster)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generierung der Netzwerkstruktur
%   im Wesentlichen wird die C2C_stat.C2C struktur genommen
%   und alle vorhandenen Konnektivitaetsmasse fuer die definierten
%   Netzwerkstrukturen umgerechnet
% output ... NET ... Struktur analog zu C2C struktur nur mit Angaben fuer
%   Netzwerke


% generiere Gruppen Verbindungsmatrix der Netzwerke GNM
NET.GNM = get_network_matrix(Cluster);
% generiere Gruppen Verbindungsmatrix zwischen den Netzwerken
[NET.GiNM,NET.inter_network_conn] = get_inter_network_matrix(Cluster);
NET.networknames = Cluster{1,1}.net.global.network_names;
NET.clusternames = Cluster{1,1}.net.global.network_names;
NET.networknumbers = [1:length(NET.networknames)];
[NET.network_indices] = get_network_names_numbers(Cluster);
[NET.GiNMq] = get_internetwork_square_GiNM(NET.GiNM,NET.inter_network_conn,length(NET.networknames));

assignin('base','GNM',NET.GNM);
assignin('base','GiNM',NET.GiNM);
assignin('base','network_indices',NET.network_indices);


% [RN] = get_network_data(C2C_stat_g1.C2C.PCC.R,NET.GNM,length(NET.networknames));
% 
% [RiN, RiNstd] = get_internetwork_data(C2C_stat_g1.C2C.PCC.R,NET.GiNM,length(NET.networknames),NET.inter_network_conn);
% 
% % 
% ZN = ck_rc3_fisher_transformation(RN);
% ZiN = ck_rc3_fisher_transformation(RiN);
% 
% NET.R = merge_net(RiN,RN);
% %C2C_stat_g1.Rstd = merge_net(RiNstd,RN);
% NET.Z = merge_net(ZiN,ZN);
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% work here further

% f = fieldnames(C2C_stat);
% for i=1:length(f)
%     if isfield(C2C_stat,f{i}) && isstruct(C2C_stat.(f{i})) 
%         
%         X.(f{i}) = adapt_C2Csubstruct_to_NETsubstruct(C2C_stat.(f{i}),NET.GNM,NET.GiNM,NET.networknames,NET.inter_network_conn,NET.network_indices);
%         
%     end 
% end

% umwandlung der C2Csubstruct in die Netzwerk substructure
Csubstruct = C2C_stat.C2C;
f = fieldnames(Csubstruct);
for i=1:length(f)
    if isfield(Csubstruct,f{i}) && isstruct(Csubstruct.(f{i}))
        f2 = fieldnames(Csubstruct.(f{i}));
        for j=1:length(f2)

            if isfield(Csubstruct.(f{i}),f2{j}) && isnumeric(Csubstruct.(f{i}).(f2{j}))
                %fprintf('entering convert_C2CMatrix_to_NETMatrix with %s.%s\n',f{i},f2{j});
              %  NETsubstruct.(f{i}) = convert_C2CMatrix_to_NETMatrix(Csubstruct.(f{i}).(f2{j}),GNM,GiNM,networknames,inter_network_conn,network_indices);
                NET.(f{i}).(f2{j}) = convert_C2CMatrix_to_NETMatrix(Csubstruct.(f{i}).(f2{j}),...
                    NET.GNM,NET.GiNM,NET.networknames,NET.inter_network_conn,NET.network_indices);
                
            end
        end
    end
end



% 
% function NETsubstruct = adapt_C2Csubstruct_to_NETsubstruct(Csubstruct,GNM,GiNM,networknames,inter_network_conn,network_indices)
% % umwandlung der C2Csubstruct in die Netzwerk substructure
% 
% f = fieldnames(Csubstruct);
% for i=1:length(f)
%     if isfield(Csubstruct,f{i}) && isstruct(Csubstruct.(f{i}))
%         f2 = fieldnames(Csubstruct.(f{i}));
%         for j=1:length(f2)
% 
%             if isfield(Csubstruct.(f{i}),f2{j}) && isnumeric(Csubstruct.(f{i}).(f2{j}))
%                 %fprintf('entering convert_C2CMatrix_to_NETMatrix with %s.%s\n',f{i},f2{j});
%               %  NETsubstruct.(f{i}) = convert_C2CMatrix_to_NETMatrix(Csubstruct.(f{i}).(f2{j}),GNM,GiNM,networknames,inter_network_conn,network_indices);
%                 NETsubstruct.(f{i}) = convert_C2CMatrix_to_NETMatrix(Csubstruct.(f{i}).(f2{j}),GNM,GiNM,networknames,inter_network_conn,network_indices);
%                 
%             end
%         end
%     end
% end


function NETR = convert_C2CMatrix_to_NETMatrix(R,GNM,GiNM,networknames,inter_network_conn,network_indices)

if length(size(R))==3
[RN] = get_network_data(R,GNM,length(networknames));
[RiN, RiNstd] = get_internetwork_data(R,GiNM,length(networknames),inter_network_conn);
NETR = merge_net(RiN,RN);
end
if length(size(R))==2
    [RN] = get_network_data2d(R,GNM,length(networknames),network_indices);
    % es gibt keine sinnvolle interpretation fuer internetwork bei 2 D
    % daten
    %[RiN, RiNstd] = get_internetwork_data2d(R,GiNM,length(networknames),inter_network_conn,network_indices);
    NETR = RN; %merge_net2d(RiN,RN);
end    


function NiN = merge_net(NiN,N)
%NiN die internetwork Matrix
% N die Network Matrix
if size(N,1)~=size(NiN,1) || length(size(N))~=2
    fprintf('error in matrix size in merge_net in ck_rc3_get_C2C_network_values2\n');
    size(N)
    size(NiN)
end
    
for i=1:size(N,1)
    NiN(i,i,:)= N(i,:);
end



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

function [RG, RGstd] = get_network_data2d(R,GM,num_nets,network_indices)

%fprintf('2d\n');
RG = zeros(num_nets,size(R,2));

RGstd = zeros(num_nets,size(R,2));

% alle Netzwerkeintraege > 0 sind dann elemente des Netzwerks

for p=1:size(R,2)
    for n=1:length(network_indices) % ueber die Netzwerke
        
        network = network_indices{n};
        
        dummyR = 0;
        for i=1:length(network) % ueber die Netzwerkmatrix
            
            dummyR = dummyR + R(network(i),p);
        end
        RG(n,p) = mean(dummyR);
        RGstd(n,p) = std(dummyR);  
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
    


