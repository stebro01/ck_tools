function estimate = ck_rc3_dentropy(xhist)
% function dentropy estimates the entropy of a time series x
% which was previously transformed in its histogram histx 
% take attention that the sum of histogram xhist is 1
% May03

[NRow,NCol]=size(xhist);

if NRow==1, xhist=transpose(xhist); end
[NRow,NCol]=size(xhist);

estimate=0;
count=0;

%NRow
for n=1:NRow
  if xhist(n)~=0 
    logterm=log2(xhist(n));
  else
    logterm=0;
  end
  count=count+xhist(n);
  estimate=estimate-xhist(n)*logterm;
end
