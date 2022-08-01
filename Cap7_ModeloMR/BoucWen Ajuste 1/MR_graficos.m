%% Gráficos del ajuste del modelo Bouc-Wen modificado a los datos experimentales
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

%% Parametros del modelo ajustado
load parametros

%% Gráficos de los resultados sin modelo 

zci = @(v) find(diff(sign(v))); % funcion para encontrar los cruces por 0

for i = 1:length(D)
    D_for = D(i);
    [fun, F_sim, F_lab, x, dx, t] = MR_BoucWen(D_for,parametros(:,i));
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
    
    subplot(1,2,2);
    plot(dx,F_lab,'DisplayName',[num2str(D_for.Amp),' A'],'linewidth',2); hold on; 
    grid on, xlabel('Velocidad [mm/s]'), ylabel('Fuerza [N]');   
    ylim([-2000 2000])
    xlim([-150 150])
    lg  = legend('Orientation','Vertical','NumColumns',1,'Location','SouthEast'); 

    if rem(i,5)==0
%         exportgraphics(gcf,strcat('Figs/',D_for.mm,"mm",D_for.freq,"Hz",'.jpg'),"Resolution",1000)
    end
end

%% Ordenar datos
amps = 0:0.25:1;
actual_plot = strcat(D(1).mm,D(1).freq);
max_forces = zeros(5,5);
errores = zeros(5,5);
param_gamma = zeros(5,5);
param_alpha = zeros(5,5);
param_c0 = zeros(5,5);
param_c1 = zeros(5,5);
count = 1;
amp_index = 1;
names = [strcat(D(1).mm, " mm, ",D(1).freq," Hz")];

for i = 1:length(D)
    D_for = D(i);
    [fun, F_sim, F_lab, x, dx, t] = MR_BoucWen(D_for,parametros(:,i));
    if ~strcmp(strcat(D(i).mm,D(i).freq), actual_plot)
        actual_plot = strcat(D(i).mm,D(i).freq);
        names = [names strcat(D(i).mm, " mm, ",D(i).freq," Hz")];
        amp_index = 1;
        count = count+1;
    end
    errores(amp_index,count) = fun;
    forces = D_for.data(:,4);
    max_forces(amp_index,count) = max(abs(forces));  
    param_alpha(amp_index,count) = parametros(1,i);
    param_gamma(amp_index,count) = parametros(2,i);
    param_c0(amp_index,count) = parametros(3,i);
    param_c1(amp_index,count) = parametros(4,i);
    amp_index = amp_index+1;
end

%% Gráfico del error RMS vs Corriente    
    
gcf = figure('Position', [10 10 600 300]);
plot(amps,errores*100,'-o','linewidth',2);
grid on, xlabel('Corriente [A]'), ylabel('Error RMS normalizado [%]');
legend(names,'Location','northeast')
xticks([0,0.25,0.5,0.75,1])
% exportgraphics(gcf,'Figs/RMSvsA1.jpg',"Resolution",1000)

%% Regresión de los parámetros
x = amps'*ones(1,size(param_gamma,2));
x = reshape(x,[],1);
y = reshape(param_gamma,[],1);

poly_gamma = fit(x,y,'exp1');
polyval_gamma = poly_gamma.a*exp(poly_gamma.b*amps);

gcf = figure('Position', [10 10 400 250]);
plot(amps,param_gamma,'-o','linewidth',1); hold on;
plot(amps,polyval_gamma,'--r','linewidth',4);
grid on, xlabel('Corriente [A]'), ylabel('Param - \gamma');
legend([names,'Ajuste'],'Location','northeast')
xticks([0,0.25,0.5,0.75,1])
title(['\gamma(i_c) = ',sprintf('%.2f',poly_gamma.a),'exp(',sprintf('%.2f',poly_gamma.b),'i_c)'])
% exportgraphics(gcf,'Figs/gamma.jpg',"Resolution",1000)

poly_alpha = polyfit(amps'*ones(1,size(param_alpha(:,[1:2,4:5]),2)),param_alpha(:,[1:2,4:5]),2);
polyval_alpha = polyval(poly_alpha,amps);

gcf = figure('Position', [10 10 400 250]);
plot(amps,param_alpha,'-o','linewidth',1); hold on;
plot(amps,polyval_alpha,'--r','linewidth',4);
grid on, xlabel('Corriente [A]'), ylabel('Param - \alpha');
% legend([names,'Ajuste'],'Location','southeast')
xticks([0,0.25,0.5,0.75,1])
title(['\alpha(i_c) = ',sprintf('%.2f',poly_alpha(1)),'i_c^2+',sprintf('%.2f',poly_alpha(2)),'i_c+',sprintf('%.2f',poly_alpha(3))])
% exportgraphics(gcf,'Figs/alpha.jpg',"Resolution",1000)

poly_c0 = polyfit(amps'*ones(1,size(param_c0,2)),param_c0,1);
polyval_c0 = polyval(poly_c0,amps);

gcf = figure('Position', [10 10 400 250]);
plot(amps,param_c0,'-o','linewidth',1); hold on;
plot(amps,polyval_c0,'--r','linewidth',4);
grid on, xlabel('Corriente [A]'), ylabel('Param - c_0');
% legend([names,'Ajuste'],'Location','northwest')
xticks([0,0.25,0.5,0.75,1])
title(['c_0(i_c) = ',sprintf('%.2f',poly_c0(1)),'i_c+',sprintf('%.2f',poly_c0(2))])
% exportgraphics(gcf,'Figs/c_0.jpg',"Resolution",1000)

poly_c1 = polyfit(amps'*ones(1,size(param_c1,2)),param_c1,1);
polyval_c1 = polyval(poly_c1,amps);

gcf = figure('Position', [10 10 400 250]);
plot(amps,param_c1,'-o','linewidth',1); hold on;
plot(amps,polyval_c1,'--r','linewidth',4);
grid on, xlabel('Corriente [A]'), ylabel('Param - c_1');
% legend([names,'Ajuste'],'Location','northwest')
xticks([0,0.25,0.5,0.75,1])
title(['c_1(i_c) = ',sprintf('%.2f',poly_c1(1)),'i_c+',sprintf('%.2f',poly_c1(2))])
% exportgraphics(gcf,'Figs/c_1.jpg',"Resolution",1000)
