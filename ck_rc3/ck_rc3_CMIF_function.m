function varargout=ck_rc3_CMIF_function(varargin)
% CMIF_function
% ursprünglich von Prof. Hoyer, dann in ck_resting_conn im 02/2011
% integriert

x=varargin{1};
y=varargin{2};
tau_range=varargin{3};
bin=varargin{4};

leng=length(x);
leng_calcu=leng-2*tau_range-1;

for tau=-tau_range:tau_range
    clear xtau; clear ytau;
    xtau=x(tau_range+1:tau_range+1+leng_calcu);
    ytau=y(tau_range+1+tau:tau_range+1+leng_calcu+tau);
    histx=hist(xtau,bin)/leng_calcu; Hx=ck_rc3_dentropy(histx);
    histy=hist(ytau,bin)/leng_calcu; Hy=ck_rc3_dentropy(histy);
    histxy=ck_rc3_dhistogram2(xtau,ytau,bin,bin)/leng_calcu;
    Hxy=ck_rc3_dentropy2(histxy);
    I = Hx + Hy - Hxy; % mutual entropy, which is the information shared by two random variables X and Y
    CMIF(tau+tau_range+1,1)=tau;
    CMIF(tau+tau_range+1,2)=I;
end
%assignin('base','histx',histx);
%assignin('base','histy',histy);
%assignin('base','histxy',histxy);
varargout{1}=CMIF;