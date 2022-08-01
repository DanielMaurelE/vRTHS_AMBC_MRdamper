function [] = sincronize(Datos,frec,Ts,amp)
% Entradas:
%     Datos : cell con todos los archivos
%     frec  : frecuencia de la sinusoidal
%     Ts    : tiempo de sampleo
%     amp   : amplitud sinusoidal

% Salidas:
%     Se guardan en la carpeta directamente los datos preprocesados

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if frec == 5 || frec == 0.5  || frec==0.1 || frec==0.25
    Data(:,2) = Datos(:,1)*980.665; % Fuerza en N
    Data(:,1) = Datos(:,2)*10; % desp en mm
else
    Data(:,2) = Datos(:,2)*980.665; % Fuerza en N
    Data(:,1) = Datos(:,1)*10; % desp en mm
end
windowSize = round(1/(10*Ts*frec));
b = (1/windowSize)*ones(1,windowSize);
a = 1;
D = filter(b,a,Data);

locs = find(isnan(Data(1:end-2,1))-isnan(Data(3:end,1)));
intercept = [400;locs;length(Data)];
cont = 1;

for i = 1:2:9
    Data_aux = Data(intercept(i):intercept(i+1),:);
    D_aux = D(intercept(i):intercept(i+1),:);
    pos = find(abs(diff(D_aux(:,1)))>Ts);
    Data_aux = D_aux(pos(1):pos(end),:);
    Data_aux = detrend(Data_aux,0);
    T = (0:Ts:(length(Data_aux)-2)*Ts)'; % 
    x = Data_aux(2:end,1);
    dx = diff(Data_aux(:,1))/Ts;
    F = Data_aux(2:end,2);
    aux = [T, x, dx, F];
    
    Amps = [0,0.25,0.5,0.75,1];
    % Escribir .txt
    fileID = [num2str(amp),'mm',strrep(num2str(frec),'.','x'),'Hz',strrep(num2str(Amps(cont)),'.','x'),'A','.txt'];
    save(['Preprocesados/',fileID], 'aux', '-ascii')
    cont = cont+1;
end

    

