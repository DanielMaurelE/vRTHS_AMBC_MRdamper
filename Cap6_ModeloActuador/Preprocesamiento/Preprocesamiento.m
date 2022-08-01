%% Datos experimentales con prensa universal uniaxial MTS sobre AMR
% Desplazamiento comandado tipo barrido sinusoidal
% Datos obtenidos por María Quiroz
% Estudiante de Magíster en Ciencias de la Ingeniería Civil, UTFSM

%% Inicializar
clear all; close all; clc
addpath('Datos');

%% Importar señales
datos = {};
% Chirp actuador solo
Data = importdata('ACT.txt').data;
Ts = 0.001; %s
datos{1} = sincronize(Data,Ts);

% Chirp 0A
Data = importdata('0A.txt').data;
Ts = 0.001; %s
datos{2} = sincronize(Data,Ts);

% Chirp con 0.5A
Data = importdata('0x5A.txt').data;
Ts = 0.001; %s
datos{3} = sincronize(Data,Ts);

% Chirp con 1A
Data = importdata('1A.txt').data;
Ts = 0.001; %s
datos{4} = sincronize(Data,Ts);

%% Guardar datos x_t,x_m y F
min_len = length(datos{1});

for i = 1:length(datos)
    if length(datos{i}) <= min_len
        min_len = length(datos{i});
        data_ref = datos{i};
    end
end

figure; 
for i = 1:length(datos)
    [i1,i2] = findsignal(datos{i}(:,1),data_ref(:,1));
    x_t(:,i) = datos{i}(i1:i2,1);
    x_m(:,i) = datos{i}(i1:i2,2);
    F_m(:,i) = datos{i}(i1:i2,3);
    plot(x_t)
end

save('frec','x_t','x_m','F_m')


