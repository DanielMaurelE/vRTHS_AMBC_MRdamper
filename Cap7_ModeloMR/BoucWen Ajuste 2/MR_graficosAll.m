%% Gráficos del ajuste del modelo Bouc-Wen modificado a los datos experimentales
% Datos obtenidos por María Quiroz
% Estudiante de Magíster en Ciencias de la Ingeniería Civil, UTFSM
clear all; close all; clc
addpath('Validacion')

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

%% Parametros del modelo ajustado
load parametros

%% Gráficos de los resultados sin modelo 

zci = @(v) find(diff(sign(v))); % funcion para encontrar los cruces por 0

for i = 1:length(D)
    D_for = D(i);
    [fun, F_sim, F_lab, x, dx, t] = MR_BoucWen(D_for,parametros);
    pos = zci(F_lab);
    index = floor(length(pos)/2);
    F_sim = F_sim(pos(index):pos(index+2));
    F_lab = F_lab(pos(index):pos(index+2));
    x = x(pos(index):pos(index+2));
    dx = dx(pos(index):pos(index+2));
    
    gcf = figure(ceil(i/5));
    gcf.Position = [10 10 800 350];
  
    subplot(1,2,1);
    plot(x,F_lab,'DisplayName',[num2str(D_for.Amp),' A'],'linewidth',2), hold on,  
    grid on, xlabel('Desplazamiento [mm]'), ylabel('Fuerza [N]'),
    ylim([-2000 2000])
%     sgtitle(strcat(D(i).mm," mm ", D(i).freq, " Hz"))
    
    subplot(1,2,2);
    plot(dx,F_lab,'DisplayName',[num2str(D_for.Amp),' A'],'linewidth',2); hold on; 
    grid on, xlabel('Velocidad [mm/s]'), ylabel('Fuerza [N]');   
    ylim([-2000 2000])
    xlim([-100 100])
    lg  = legend('Orientation','Vertical','NumColumns',1,'Location','SouthEast'); 

    if rem(i,5)==0
        exportgraphics(gcf,strcat('Figs/',D_for.mm,"mm",D_for.freq,"Hz",'.jpg'),"Resolution",1000)
    end
end

%% Ordenar datos
amps = 0:0.25:1;
actual_plot = strcat(D(1).mm,D(1).freq);
max_forces = zeros(5,5);
errores = zeros(5,5);
count = 1;
amp_index = 1;
names = [strcat(D(1).mm, " mm, ",D(1).freq," Hz")];

for i = 1:length(D)
    D_for = D(i);
    [fun, F_sim, F_lab, x, dx, t] = MR_BoucWen(D_for,parametros);
    if ~strcmp(strcat(D(i).mm,D(i).freq), actual_plot)
        actual_plot = strcat(D(i).mm,D(i).freq);
        names = [names strcat(D(i).mm, " mm, ",D(i).freq," Hz")];
        amp_index = 1;
        count = count+1;
    end
    errores(amp_index,count) = fun;
    forces = D_for.data(:,4);
    max_forces(amp_index,count) = max(abs(forces));  
    amp_index = amp_index+1;
end
        
%% Gráfico de la Fuerza vs Corriente    
Markers = '*+vxd';
gcf = figure('Position', [10 10 600 300]);
for i = 1:5
    midx = 5-rem(i,length(Markers));   %cycle through them
    plot(amps,max_forces(:,i),':','Marker', Markers(midx),'linewidth',2);
    hold on;
end
grid on, xlabel('Corriente [A]'), ylabel('Fuerza [N]');
legend(names,'Location','southeast')
xticks([0,0.25,0.5,0.75,1])
exportgraphics(gcf,'Figs/FvsA.jpg',"Resolution",1000)

%% Gráfico del error RMS vs Corriente      
gcf = figure('Position', [10 10 600 300]);
plot(amps,errores*100,'-o','linewidth',2);
grid on, xlabel('Corriente [A]'), ylabel('Error RMS normalizado [%]');
legend(names,'Location','southeast','NumColumns',3)
xticks([0,0.25,0.5,0.75,1])
ylim([0 18])
exportgraphics(gcf,'Figs/RMSvsA.jpg',"Resolution",1000)

%% Graficos con modelo ajustado
zci = @(v) find(diff(sign(v))); % funcion para encontrar los cruces por 0

for i = 1:length(D)  
    [fun, F_sim, F_lab, x, dx, t] = MR_BoucWen(D(i),parametros); 
    pos = zci(F_lab);
    index = floor(length(pos)/2);
    F_sim = F_sim(pos(index):pos(index+2));
    F_lab = F_lab(pos(index):pos(index+2));
    x = x(pos(index):pos(index+2));
    dx = dx(pos(index):pos(index+2));
    
    gcf = figure(ceil(i/5)+11);
    gcf.Position = [10 10 800 350];

    subplot(1,2,1);
    plot(x,F_lab,'k','linewidth',2), hold on,
    plot(x,F_sim,'--r','linewidth',2), 
    grid on, xlabel('Desplazamiento [mm]'), ylabel('Fuerza [N]'),
    ylim([-2000 2000])
    
    subplot(1,2,2);
    plot(dx,F_lab,'k','linewidth',2); hold on; 
    plot(dx,F_sim,'--r','linewidth',2),

    grid on, xlabel('Velocidad [mm/s]'), ylabel('Fuerza [N]');
    xlim([-150 150])
    ylim([-2000 2000])
    lg  = legend('Experimental','Simulación');
    lg.Location = 'SouthEast'; 
%     sgtitle(strcat(D(i).mm," mm ", D(i).freq, " Hz"))
    if rem(i,5)==0
        exportgraphics(gcf,strcat('Figs/',D(i).mm,"mm",D(i).freq,"Hz",'-modelo','.jpg'),"Resolution",1000)
    end
end

%% Importar datos de validación
val = load('Validacion/validacion').data;
VAL = struct();
cor_val = [0.2:0.2:1];

for i = 1:size(val,2)-1
    Ts = 0.001;
    x_val = val(2:end,1);
    dx_val = diff(val(:,1))/Ts;
    f_val = val(2:end,i+1);
    t_val = [0:0.001:(length(val)-1)*0.001]';
    t_val = t_val(1:end-1);
    VAL(i).data = [t_val,x_val,dx_val,f_val];
    VAL(i).Amp = cor_val(i);
end

%% Calcular error datos validacion
cont = 2;
val_sim = [];
val_error = zeros(1,size(VAL,2));
for i = 1:size(VAL,2)
    [fun, F_sim, F_lab, x, dx, t] = MR_BoucWen(VAL(i),parametros);    
    val_sim(:,i) = F_sim;
    val_error(1,i) = fun;
end

%% Graficos validación
% Desplazamiento
gcf = figure('Position', [10 10 900 200]);
plot(t_val,x_val,'k','linewidth',1);
grid on, xlabel('Tiempo [s]'), ylabel('Desplazamiento [mm]');
xlim([0,t_val(end)])
exportgraphics(gcf,'Figs/desp_val.jpg',"Resolution",1000)

% Fuerzas
limy = [400,800,1000,1500,1500];
for i = 1:size(VAL,2)
    gcf = figure('Position', [10 10 900 200]);
    plot(VAL(i).data(:,1),VAL(i).data(:,4),'k','linewidth',1,'DisplayName','Experimental'); hold on;
    plot(t,val_sim(:,i),'--r','linewidth',1,'DisplayName','Simulación');
    grid on, xlabel('Tiempo [s]'), ylabel('Fuerza [N]');
    xlim([0,t_val(end)])
    ylim([-limy(i) limy(i)])
    legend('Location','southeast')
    exportgraphics(gcf,['Figs/fuerza_val',num2str(i),'.jpg'],"Resolution",1000)
end

% Error
gcf = figure('Position',[10 10 600 300]);
plot(cor_val,val_error.*100,'o-k','linewidth',1)
xticks(0:0.2:1)
% xlim([0 1])
grid on;
xlabel('Corriente [A]'), ylabel('Error RMS normalizado [%]'),
exportgraphics(gcf,'Figs/Error_validacion.jpg',"Resolution",1000)



