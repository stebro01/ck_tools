function C = ck_rc3_estimate_C2C_statistic2(outdir,C1,C2, numgroup, method, normalize)
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
    C1 = normalize_intersubject_Z_values(C1);
    C2 = normalize_intersubject_Z_values(C2);
end

C = struct; % C2C_struct Struktur
if  numgroup==2
    f = fieldnames(C1);
    for i=1:length(f)
        if isfield(C1,f{i}) && isfield(C2,f{i}) && isstruct(C1.(f{i})) && isstruct(C2.(f{i}))
            %fprintf('entering estimate_C2C with %s %s\n',f{i},f{i});
            C.(f{i}) = estimate_C2C(method,C1.(f{i}),C2.(f{i}));
            
        elseif isfield(C1,f{i}) && isfield(C2,f{i})
            
            C.([f{i} '1'])=C1.(f{i});
            C.([f{i} '2'])=C2.(f{i});
            
        end
        
        
    end
    
end



function C = estimate_C2C(method,C1,C2)
% vereinigt die C2C_stat.C2C strukturen von beiden Gruppen
C = struct;
f = fieldnames(C1);
for i=1:length(f)
    
    if isfield(C1,f{i}) && isfield(C2,f{i}) && isstruct(C1.(f{i})) && isstruct(C2.(f{i}))  
        fprintf('entering merge_conn_field with %s %s\n',f{i},f{i});
        C.(f{i}) = merge_conn_field(C1.(f{i}),C2.(f{i}));
    elseif isfield(C1,f{i}) && isfield(C2,f{i}) 
         C.(f{i}) = C1.(f{i});
    end
    
end    





function PCC = merge_conn_field(PCC1,PCC2)
% vereinigt die PCC strukturen
PCC = struct;
%fprintf('fuerher nuen folgede strukturen zusammen\n');
f = fieldnames(PCC1);
for i=1:length(f)
    PCC = mergex(PCC,PCC1,PCC2,f{i});
end    
assignin('base','PCC',PCC);
% PCC = mergex(PCC,PCC1,PCC2,'R');
% PCC = mergex(PCC,PCC1,PCC2,'Z');
% PCC = mergex(PCC,PCC1,PCC2,'P');


function C = mergex(C,C1,C2,str)

if isfield(C1,str) && isfield(C2,str) && isnumeric(C1.(str)) && isnumeric(C2.(str))
    C.([str '1']) = C1.(str);
    C.([str '2']) = C2.(str);
    strnew = ['MeanGroupDiff' str];
    C.(strnew) = estimate_MeanGroupDiff(C1.(str),C2.(str));
end
    



function MeanGroupDiff = estimate_MeanGroupDiff(X1,X2)
if length(size(X1))==3 && length(size(X2))==3
    for i=1:size(X1,1)
        for j=1:size(X1,2)
            if i==j
                MeanGroupDiff(i,j)=0;
            else
                MeanGroupDiff(i,j) = mean(squeeze(X1(i,j,:)))-mean(squeeze(X2(i,j,:)));
            end
        end
    end
elseif length(size(X1))==2 && length(size(X2))==2
    for i=1:size(X1,1)
         MeanGroupDiff(i) = mean(squeeze(X1(i,:)))-mean(squeeze(X2(i,:)));
    end
end



function C = normalize_intersubject_Z_values(C)
% Die Datensätze werden um den Absolutbetrag der summe aller R-Werte
% korrigiert um z.B. altersabhängige Steigerungen oder Senkungen der
% Connektivität die lokal unspezifisch sind zu kompensieren
%  dann könnte ich aber keine R-Z-Transformation mehr durchführen
% also korrigiere ich einfach nur die Z-Werte
%data=data;
% ermittle den durchschnittlichen Wert der Z-Werte
f = fields(C);
for i=1:length(f)
    if isstruct(C.(f{i}))
        f2 = fields(C.(f{i}));
        for j = 1:length(f2)
            if isstruct(C.(f{i}).(f2{j}))
                if isfield(C.(f{i}).(f2{j}),'Z')
                    C.(f{i}).(f2{j}).Z = normalize_Z(C.(f{i}).(f2{j}).Z);
                end
            end
        end
    end
end




function Z = normalize_Z(Z)

%%%%%SB ADDED %%%%%%%%%%%%
%remove infinities >> aus irgendeinem Grund stehen in den Diagonalen INF
%und sehr große Werte >> die werden genullt, dann läuft die
%Normalize-Routine
for i=1:size(Z,3)
    for j = 1:size(Z,2)
        Z(j,j,i) = 0;
        Z(Z == Inf) = 0;
    end
end
%%%%END

for i=1:size(Z,3)
    x(i)=sum(sum(abs(Z(:,:,i))));
end
mr=mean(x)/(size(Z,1)*(size(Z,2)-1));
m=mean(x);

for i=1:size(Z,3)
    X1=squeeze(Z(:,:,i));
    mi=sum(sum(abs(X1)));
    X1=X1*(m/mi);
    Z(:,:,i)=X1;
end

%
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


