function aux = sincronize(Datos,Ts)
% Entradas:
%     Datos : cell con todos los archivos
%     Ts    : tiempo de sampleo

% Salidas:
%     x     : desplazamiento medido
%     F     : Fuerza medida
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Data(:,2) = Datos(:,2)*980.665; % Fuerza en N
Data(:,1) = Datos(:,1)*10; % desp en mm

[~,locs] = findpeaks(Data(:,1),'MinPeakProminence',1);

windowSize = 300; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
D = filter(b,a,Data);
D = D(locs(1)-1000:end,:);
D = D-D(end,:);
D(:,2) = detrend(D(:,2),0);

x = D(:,1);
F = D(:,2);
aux = [x, F];