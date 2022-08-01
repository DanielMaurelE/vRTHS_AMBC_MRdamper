%% Ajuste del modelo Bouc-Wen modificado a los datos experimentales
% Datos obtenidos por María Quiroz
% Estudiante de Magíster en Ciencias de la Ingeniería Civil, UTFSM
clear all; close all; clc

%% Data de entrenamiento
myDir = 'Preprocesamiento-sinusoidales/Preprocesados/';
myFiles = dir(fullfile(myDir,'*.txt'));
D = struct();
File = 1;
for i = 1:length(myFiles)
    data_name = split(myFiles(i).name,'mm');
    mm = data_name{1};
    
    data_name = split(data_name{2},["Hz","A.txt"]);
    Ampere = strrep(data_name{2},'x','.');
    freq = strrep(data_name{1},'x','.');
    
    D(File).mm = mm;
    D(File).freq = freq;
    D(File).Amp = str2double(Ampere);
    D(File).data = importdata([myDir,myFiles(i).name]);
    File = File+1;    
end


%% Optimizacion
fun = @(params) MR_BoucWenAll(D,params);

lb = [0,0,0,-inf*ones(1,9)];
ub = [];
x0 = [0,0,0,-17.1,23.82,5.73,0.68,-2.91,12.34,0.67,86.19,10.4];

options = optimoptions('fmincon','Display','iter');
param = fmincon(fun,x0,[],[],[],[],lb,ub,[],options);
parametros = param';
save('parametros','parametros');

for i = 1:length(D)
    [fun, F_sim, F_lab, x, dx, t] = MR_BoucWen(D(i),parametros);
    fun
    gcf = figure(i);
    set(gcf,'position',[10,10,1000,600])
    subplot(2,2,[1 2]), plot(t,F_lab,'k','DisplayName',[num2str(D(i).Amp),' A']), hold on, plot(t,F_sim,'--r'), grid on, xlabel('Tiempo [s]'), ylabel('F [N]'), xlim([0 t(end)]);
    subplot(2,2,3), plot(x,F_lab,'k','DisplayName',[num2str(D(i).Amp),' A']), hold on, plot(x,F_sim,'--r'), grid on, xlabel('Desp [mm]'), ylabel('F [N]'), 
    subplot(2,2,4), plot(dx,F_lab,'k','DisplayName',[num2str(D(i).Amp),' A']), hold on, plot(dx,F_sim,'--r'), grid on, xlabel('Velocidad [mm/s]'), ylabel('F [N]');
    sgtitle(strcat(D(i).mm," mm ", D(i).freq, " Hz"))
end







