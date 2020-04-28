
% % Multiregressionsanalyse
% function [reg] = estimate_multiregression(Tab,titel)
target = 'door2needle';
col_doo =  find(strcmp(titel,target));

col_tra = find(strcmp(titel,'Transp'));
col_med = find(strcmp(titel,'MedRD'));
col_ver = find(strcmp(titel,'VerdRDgrupp'));
col_nii = find(strcmp(titel,'NIHSSin'));
col_ent=find(strcmp(titel,'EntfernungRettungsweg'));

% hole reine Zahlenspalten aus den spalten
xdoo = convert_vector_to_numbers(Tab(:,col_doo),4,270);
assignin('base','xdoo',xdoo);
xtra = convert_vector_to_numbers(Tab(:,col_tra),0,10);
assignin('base','xtra',xtra);
xmed = convert_vector_to_numbers(Tab(:,col_med),0,1);
assignin('base','xmed',xmed);
xver = convert_vector_to_numbers(Tab(:,col_ver),0,20);
assignin('base','xver',xver);
xnii = convert_vector_to_numbers(Tab(:,col_nii),0,50);
assignin('base','xnii',xnii);
xent = convert_vector_to_numbers(Tab(:,col_ent),-1,70);
assignin('base','xent',xent);

% baue nun die Matrix zusammen
% korrigiere noch unabhaengigen Variablen 
% und wandele sie in Dummy variablen um siehe S 13 in Multivariate
% Analysemethoden

% Transport 0 RD 1 NA sonst NaN
xtra(xtra==0)=NaN;
xtra(xtra==1)=0; % RD
xtra(xtra==2)=1;
xtra(xtra==3)=1;
xtra(xtra>3)=NaN;

xmed(xmed>1)=NaN;
%Verdachtsdiagnose  1.. Stroke 0 sonst
xver(xver>1)=0;

% Fuege Matrix zusammen
y = xdoo';
X = xtra';
X(:,2) = xmed;
X(:,3) = xver;
X(:,4) = xnii;
X(:,5) = xent;
X(:,6) = 1;

assignin('base','y',y);
assignin('base','X',X);

% b koeffizienten
% bint a p-by-2 matrix bint of 95% confidence intervals for the coefficient estimates
% r residuals
% rint   interval that can be used to identify outliers
% stats ... a 1-by-4 vector stats.. R2 statistic, the F statistic and its p value, and an estimate of the error variance.
%[b,bint,r,rint,stats] = regress(y,X);
lm = fitlm(X,y,'linear')
reg = 0;

end
