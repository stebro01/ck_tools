function estimate = ck_rc3_dentropy2(xyhist)
% function dentropy2 estimates the COMMON entropy of two time series x,y
% which was previously transformed in its xy-histogram xyhist
% take attention that the sum of histogram xyhist is 1
% May03

[NRow,NCol]=size(xyhist);

estimate=0;
count=0;

for n=1:NCol
    for m=1:NRow
        if xyhist(n,m)~=0 
        logterm=log2(xyhist(n,m));
        else
        logterm=0;
        end;
    count=count+xyhist(n,m);
    estimate=estimate-xyhist(n,m)*logterm;
    end;
end;

estimate=estimate/count;
