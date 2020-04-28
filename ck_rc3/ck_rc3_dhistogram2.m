function [result]=ck_rc3_dhistogram2(x,y,ncellx,ncelly)
% DHISTOGRAM2_may Computes the two dimensional histogram of two raw vectors x and y.
% May03

[NRowX,NColX]=size(x);
if NRowX~=1
  x=x';
end;

[NRowY,NColY]=size(y);
if NRowY~=1
  y=y';
end;

[NRowX,NColX]=size(x);
[NRowY,NColY]=size(y);

if NColX~=NColY
  error('Unequal length of X and Y');
end;

lowerx=min(x)- (max(x)-min(x))/(length(x)-1)/2; 
upperx=max(x)+ (max(x)-min(x))/(length(x)-1)/2; 
lowery=min(y)- (max(y)-min(y))/(length(y)-1)/2; 
uppery=max(y)+ (max(y)-min(y))/(length(y)-1)/2;


result(1:ncellx,1:ncelly)=0;

xx=round( (x-lowerx)/(upperx-lowerx)*ncellx + 1/2 );
yy=round( (y-lowery)/(uppery-lowery)*ncelly + 1/2 );
for n=1:NColX
  indexx=xx(n);
  indexy=yy(n);
  if indexx >= 1 & indexx <= ncellx & indexy >= 1 & indexy <= ncelly
    result(indexx,indexy)=result(indexx,indexy)+1;
  end;
end;