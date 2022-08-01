%% Preprocesamiento de los datos de validacion
% Desplazamiento comandado tipo terremoto
% Datos obtenidos por María Quiroz
% Estudiante de Magíster en Ciencias de la Ingeniería Civil, UTFSM

%% Inicializar
clear all; close all; clc
addpath('MRDAMPER')

%% Importar señales
datos = {};
Data = importdata('terremoto0x2A.txt').data;
Ts = 0.001; %s
datos{1} = sincronize(Data,Ts);
s0 = datos{1};

Data = importdata('terremoto0x4A.txt').data;
Ts = 0.001; %s
datos{2} = sincronize(Data,Ts);

Data = importdata('terremoto0x6A.txt').data;
Ts = 0.001; %s
datos{3} = sincronize(Data,Ts);

Data = importdata('terremoto0x8A.txt').data(69936:end,:);
Ts = 0.001; %s
datos{4} = sincronize(Data,Ts);

Data = importdata('terremoto1A.txt').data;
Ts = 0.001; %s
datos{5} = sincronize(Data,Ts);

%% Guardar datos [desp,fuerzas]
min_len = length(datos{1});

for i = 1:length(datos)
    if length(datos{i}) <= min_len
        min_len = length(datos{i});
        data = zeros(length(datos{i}),length(datos)+1);
        data(:,[1,i+1]) = datos{i};
        data_ref = datos{i};
    end
end


for i = 1:length(datos)
    [i1,i2] = findsignal(datos{i}(:,1),data_ref(:,1));
    data(:,i+1) = datos{i}(i1:i2,2);
    figure; plot(datos{i})
end

save('validacion','data')


