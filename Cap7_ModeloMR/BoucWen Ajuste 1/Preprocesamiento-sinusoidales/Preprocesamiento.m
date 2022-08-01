%% Datos experimentales con prensa universal uniaxial MTS sobre AMR
% Datos obtenidos por María Quiroz
% Estudiante de Magíster en Ciencias de la Ingeniería Civil, UTFSM
clear all; close all; clc
addpath('Datos_crudos');

%% Guardar datos pre-procesados
myDir = 'Datos_crudos/'; %gets directory
myFiles = dir(fullfile(myDir,'*.txt'));

for i = 1:length(myFiles)
    nombre = myFiles(i).name;
    div = split(nombre,'mm');
    k = regexp(div(1),'\d*','Match');
    amp = str2double(k{1});
    m = regexp(div(2),'\d*','Match');
    if length(m{1}) == 1
        frec = str2double(m{1});
    else
        frec = str2double([m{1}{1},'.',m{1}{2}]);
    end
    if frec == 5 || frec == 0.5 || frec==0.1 || frec==0.25
        Ts = 0.005;
    else
        Ts = 0.001;
    end
    Data = importdata(myFiles(i).name).data;
    sincronize(Data,frec,Ts,amp)
    
end

%% Figuras de los datos pre-procesados
myDir = 'Preprocesados/'; %gets directory
myFiles = dir(fullfile(myDir,'*.txt'));
D = {};
File = 1;
for i = 1:length(myFiles)
    D{File} = importdata([myDir,myFiles(i).name]);
    % Graficos
    figure()
    subplot(2,2,[1 2]), plot(D{File}(:,1),D{File}(:,4)), grid on, xlabel('Tiempo [s]'), ylabel('F [N]');
    subplot(2,2,3), plot(D{File}(:,2),D{File}(:,4)), grid on, xlabel('Desp [mm]'), ylabel('F [N]'); 
    subplot(2,2,4), plot(D{File}(:,3),D{File}(:,4)),grid on, xlabel('Velocidad [mm/s]'), ylabel('F [N]');
    sgtitle(myFiles(i).name);
    File = File+1;    
end





