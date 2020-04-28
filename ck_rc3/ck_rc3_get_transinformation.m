function CMIF = ck_rc3_get_transinformation(x,y,tau,bin)
% Berechnet die Transinformation
% x ist ein mx1 Vector,
% Y ist ein mxn Vector n... Clusteranzahl

x=detrend(x);
y=detrend(y);


%CMIF
%clear AMIFx;
AMIFx=ck_rc3_CMIF_function(x,x,tau,bin);
%plot(AMIFx(:,1),AMIFx(:,2));ylabel('AMIFx')
max_x=max(AMIFx(:,2));
%clear AMIFy
AMIFy=ck_rc3_CMIF_function(y,y,tau,bin);
%plot(AMIFy(:,1),AMIFy(:,2));ylabel('AMIFy')
max_y=max(AMIFy(:,2));
normmax=min(max_x,max_y);


% CMIF=CMIF_function(x,y,tau,bin);
% subplot(4,2,3)
% plot(CMIF(:,1),CMIF(:,2));ylabel('CMIF')

% normiert auf Gesamtinformation

CMIF=ck_rc3_CMIF_function(x,y,tau,bin);
size(CMIF);
CMIF(:,2)=CMIF(:,2)./normmax;

%subplot(4,2,3)
%plot(CMIF(:,1),CMIF(:,2));ylabel('CMIF');axis([-tau tau 0 1]);
