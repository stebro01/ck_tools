function [Cxy_bandpass F_bandpass] = ck_rc3_get_coherence(x,y,TR,f1,f2)
% berechung des Coherence zwischen x und y im Frequenzbereich 
% zwischen f1 und f2
% Es wird dabei angenommen, dass sich der Signalverlauf in Zeilenrichtung
% bewegt.
% x.... n x 1 Matrix (n Anzahl der Messungen, m Anzahl der zu
%             korrelierenden Signalverläufe)
% y.... n x 1 Matrix (n Anzahl der Messungen)
%


% berechne den Coherence
% Überprüfe ob die Dimensionsanforderungen eingehalten wurden
if  size(size(y),2)~=2 || min(size(y))~=1
    fprintf('Fehler in ck_rc_get_coherence ... die Dimension der INputvariable y muss nx1 sein');
end    
if  size(size(x),2)~=2 || min(size(x))~=1
    fprintf('Fehler in ck_rc_get_coherence ... die Dimension der INputvariable x muss nx1 sein');
end    
if size(y,2)>size(y,1)
    y=y';
end
if size(x,2)>size(x,1)
    x=x';
end




x=detrend(x);
y=detrend(y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% berechne die Coherence
% Parameterermittlung zur Coherenceberechnung
% Window
% ein längeres Window ermöglicht eine genauere FFT
% und damit genauere Frequenzaufschlüsselungen jedoch
% um den Preis einer geringeren Reliabilität
p = size(y,1);
if p>300
    window = 200;
else
    window = p-50;
end
noverlap = 90; % Overlap der windowsegmente in Prozent
nfft = 500; % Anzahl der berechneten Ergebnispunkte = nfft/2(+1)
%fprintf('size(X)=%d,%d   --- size(Y)=%d,%d \n',size(X,1),size(X,2),size(Y,1),size(Y,2));
%[Cxy,F] = mscohere(x,y,window,noverlap,500,1/TR);
[Cxy,F] = mscohere(x,y,[],[],[],1/TR);
% Schränke die Frequenz auf das gewünschte Maß ein
Cxy_bandpass = Cxy(F(F<f2)>f1);
F_bandpass = F(F(F<f2)>f1);
%for j=1:size(Y,1)
%         R = zeros(size(Y,1),size(X,2));
%         P = zeros(size(Y,1),size(X,2));
%         [Rtmp, Ptmp] = corrcoef([X Y']);
%         R = zeros(size(Y,1),1);
%         P = zeros(size(Y,1),1);
%         R(:,1) = Rtmp(2:end,1);
%         P(:,1) = Ptmp(2:end,1);
%         for i=1:size(R,1)
%             if isnan(R(i))
%                 R(i) = 0;
%                 P(i) = 1;
%             end
%         end


