function [R P] = ck_rc3_get_pearsons_correlations(X,Y)
% berechung des Pearson Korrelationskoeffizienten zwischen X und Y
% Es wird dabei angenommen, dass sich der Signalverlauf in Zeilenrichtung
% bewegt.
% X.... 1 x n Matrix (n Anzahl der Messungen, m Anzahl der zu
%             korrelierenden Signalverläufe)
% Y.... beliebige Dimension (2-4) die letze Dimension muss dabei n betragen
%       2-Dim p x n
%       3-Dim p x q x n
%       4-Dim p x q x r x n
%
% R.... p x q x r x 1 is the zeroth lag of the normalized covariance function
% P.... p x q x r x 1 matrix of p-values for testing the hypothesis of
%       no correlation. Each p-value is the probability of getting
%       a correlation as large as the observed value by random chance,
%       when the true correlation is zero. If P(i,j) is small,
%       say less than 0.05, then the correlation R(i,j) is significant.

% berechne den Korrelationskoeffizenten
dimX = size(X);
if dimX(1)>1 && dimX(2)==1
    X = X';
end
assignin('base','X',X)
assignin('base','Y',Y)

dim = size(size(Y),2);

%assignin('base','Y',Y);
switch dim
    % Y = p x n
    case 2
        % Prüfe ob die Dimensionsanforderungen eingehalten werden
        nx = size(X,2);
        p = size(Y,1);
        ny = size(Y,2);
        if nx~=ny
            disp('Matrix size dont agree');
            if nx==p
                disp('Transpose Matrix Y');
                Y = Y';
            else
                msgbox('Error: Matrixsize dont agree');
                %PC{i,1} = -1;
                return
            end
        end
        % berechne die Korrelation
        %for j=1:size(Y,1)
        %         R = zeros(size(Y,1),size(X,2));
        %         P = zeros(size(Y,1),size(X,2));
        [Rtmp, Ptmp] = corrcoef([X' Y']);
        R = zeros(size(Y,1),1);
        P = zeros(size(Y,1),1);
        R(:,1) = Rtmp(2:end,1);
        P(:,1) = Ptmp(2:end,1);
        for i=1:size(R,1)
            if isnan(R(i))
                R(i) = 0;
                P(i) = 1;
            end
        end

    case 3
        % Prüfe ob die Dimensionsanforderungen eingehalten werden
        nx = size(X,2);
        p = size(Y,1);
        q = size(Y,2);
        ny = size(Y,3);
        if nx~=ny
            disp('Matrix size dont agree');
            msgbox('Error: Matrixsize dont agree');
            %PC{i,1} = -1;
            return
        end
        % berechne die Korrelation
        R = zeros(size(Y,1),size(Y,2));
        P = zeros(size(Y,1),size(Y,2));
        for j=1:size(Y,2)
            tmp = squeeze(Y(:,j,:));
            [Rtmp Ptmp] = corrcoef([X tmp']);
            R(:,j) = Rtmp(2:end,1);
            P(:,j) = Ptmp(2:end,1);
        end
        for i=1:size(R,1)
            for j=1:size(R,2)
                if isnan(R(i,j))
                    R(i,j) = 0;
                    P(i,j) = 1;
                end
            end
        end

    case 4
        % Prüfe ob die Dimensionsanforderungen eingehalten werden
        nx = size(X,2);
        p = size(Y,1);
        q = size(Y,2);
        r = size(Y,3);
        ny = size(Y,4);
        if nx~=ny
            disp('Matrix size dont agree');
            msgbox('Error: Matrixsize dont agree');
            %PC{i,1} = -1;
            return
        end
        % berechne die Korrelation
        R = zeros(size(Y,1),size(Y,2),size(Y,3));
        P = zeros(size(Y,1),size(Y,2),size(Y,3));
        %h = waitbar(0,'please wait');
        for j=1:size(Y,2)
            %waitbar(j/size(Y,2));
            for k=1:size(Y,3)
                tmp = squeeze(Y(:,j,k,:));
                [Rtmp Ptmp] = corrcoef([X tmp']);

                R(:,j,k) = Rtmp(2:end,1);
                P(:,j,k) = Ptmp(2:end,1);
            end
        end
        %close(h);
        for i=1:size(R,1)
            for j=1:size(R,2)
                for k=1:size(R,3)
                    if isnan(R(i,j,k))
                        R(i,j,k) = 0;
                        P(i,j,k) = 1;
                    end
                end
            end
        end
end
% Elemente mit NaN (z.B. weil alle Elemente von Y =0 sind (auserhalb des Hirns))
% werden auf 0 gesetzt für keine Korrelation -> P=1


