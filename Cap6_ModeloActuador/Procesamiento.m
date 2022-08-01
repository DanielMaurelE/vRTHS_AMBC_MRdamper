%% Análisis de los datos experimentales de prensa uniaxial MTS sobre MR
% Desplazamiento comandado tipo barrido sinusoidal
% Datos obtenidos por María Quiroz
% Estudiante de Magíster en Ciencias de la Ingeniería Civil, UTFSM

%% Inicializar
clear all; close all; clc
addpath('Preprocesamiento');
load('frec')

%% FRFs desplazamiento del actuador
% x_m vs x_t
nfft = 2^14;
w = [];
fs = 1000;

[H_1,coh_1,fuu] = H1est(x_m(:,1),x_t(:,1),w,nfft,fs);
[H_2,coh_2,~] = H1est(x_m(:,2),x_t(:,2),w,nfft,fs);
[H_3,coh_3,~] = H1est(x_m(:,3),x_t(:,3),w,nfft,fs);
[H_4,coh_4,~] = H1est(x_m(:,4),x_t(:,4),w,nfft,fs);

%% Figuras FRFs
grayColor = [.7 .7 .7];
gcf = figure('Position', [10 10 800 400]);

subplot(3,1,1);
plot(fuu,20*log10(abs(H_1)),'k','Linewidth',2,'DisplayName','Actuador solo'); hold on; 
plot(fuu,20*log10(abs(H_2)),'--', 'Color', grayColor,'Linewidth',2,'DisplayName','0 A'); 
% plot(fuu,20*log10(abs(H_3)),'r--','Linewidth',2,'DisplayName','0.5 A'); 
plot(fuu,20*log10(abs(H_4)),'r--','Linewidth',2,'DisplayName','1 A'); 
xlim([0 20]);
ylabel('Magnitud [dB]'); 
legend();
% title('Diagrama de Bode')
grid on

subplot(3,1,2);
plot(fuu,wrapTo180(angle(H_1).*-180/pi),'k','Linewidth',2,'DisplayName','Actuador solo'); hold on; 
plot(fuu,wrapTo180(angle(H_2).*-180/pi),'--','Color', grayColor,'Linewidth',2,'DisplayName','0 A'); 
% plot(fuu,wrapTo180(angle(H_3).*-180/pi),'r--','Linewidth',2,'DisplayName','0.5 A');
plot(fuu,wrapTo180(angle(H_4).*-180/pi),'r--','Linewidth',2,'DisplayName','1 A');
xlim([0 20]);
% ylim([-180,180]); 
% yticks(-180:45:180); 
ylabel('Fase [°]');  
% legend();
grid on

subplot(3,1,3);
plot(fuu,coh_1,'k','Linewidth',2,'DisplayName','Actuador solo'); hold on; 
plot(fuu,coh_2,'--','Color', grayColor,'Linewidth',2,'DisplayName','0 A'); 
% plot(fuu,coh_3,'r--','Linewidth',2,'DisplayName','0.5 A');
plot(fuu,coh_4,'r--','Linewidth',2,'DisplayName','1 A');

xlim([0 20]);
ylim([0.8 1.1]); 
% yticks(-180:45:180); 
ylabel('Coherencia');
xlabel('Frecuencia [Hz]');
% legend();
grid on
% sgtitle('x_m vs x_t')
% exportgraphics(gcf,'Figs/Frec.jpg',"Resolution",1000)


% Zoom en el retraso
tiempo = 0:1/1000:(length(x_m)-1)/1000;

gcf = figure('Position', [10 10 400 200]);

plot(tiempo,x_t(:,1),'k','Linewidth',1,'DisplayName','x_t'); hold on; 
plot(tiempo,x_m(:,1),'r--','Linewidth',1,'DisplayName','x_m'); 
xlim([4.668 6.553]); 
ylabel('Desplazamiento [mm]');
xlabel('Tiempo [s]');
legend('Location','NorthEast');
grid on
% exportgraphics(gcf,'Figs/DvsT.jpg',"Resolution",1000)


%% Identificación de sistemas
%% importar MFDID
addpath('MFDID')
mfdid

%% Modelo orden 2
load modelo
MFDID = sys_mfdid;
MFDID = idtf(MFDID);

[mag,phase,wout0] = bode(MFDID,[0.1:0.1:20*2*pi]);
mag = squeeze(mag);
phase = squeeze(phase);
wout0 = wout0./2/pi;

gcf = figure('Position', [10 10 800 400]);
subplot(3,1,1);
plot(fuu,db(abs(H_1)),'k','Linewidth',2,'DisplayName','Datos experimentales'); hold on;
plot(wout0,db(abs(mag)),'r--','Linewidth',2,'DisplayName','Modelo')
xlim([0 20])
grid on;
ylabel('Magnitud [dB]');
legend();

subplot(3,1,2);
plot(fuu,wrapTo180(angle(H_1).*-180/pi),'k','Linewidth',2,'DisplayName','Datos experimentales'); hold on;
plot(wout0,wrapTo180(phase),'r--','Linewidth',2,'DisplayName','Modelo');
xlim([0 20])
grid on;
ylabel('Fase [°]');

subplot(3,1,3);
plot(fuu,angle(H_1).*-180/pi/360./fuu.*1000,'k','Linewidth',2,'DisplayName','Datos experimentales'); hold on;
plot(wout0,phase./360./wout0.*1000,'r--','Linewidth',2,'DisplayName','Modelo');
xlim([0 20])
ylim([-50 0])
grid on;
ylabel('Retraso [ms]');
xlabel('Frecuencia [Hz]');
% exportgraphics(gcf,'Figs/sysID_frec.jpg',"Resolution",1000)

%% Modelo solo

gcf = figure('Position', [10 10 800 400]);
subplot(3,1,1);
plot(wout0,db(abs(mag)),'r','Linewidth',2); 
xlim([0 20])
grid on;
ylabel('Magnitud [dB]');

subplot(3,1,2);
plot(wout0,wrapTo180(phase),'r','Linewidth',2');
xlim([0 20])
grid on;
ylabel('Fase [°]');

subplot(3,1,3);
plot(wout0,phase./360./wout0.*1000,'r','Linewidth',2);
xlim([0 20])
ylim([-50 0])
grid on;
ylabel('Retraso [ms]');
xlabel('Frecuencia [Hz]');
% exportgraphics(gcf,'Figs/act_initialmodel.jpg',"Resolution",1000)

%% Gráfico de polos
pol = pole(MFDID);

gcf = figure('Position', [10 10 400 300]);
plot(real(pol(:)),imag(pol(:)),'rx','MarkerSize',20,'LineWidth', 2)
sgrid
xlabel('Eje real');
ylabel('Eje imaginario');
% exportgraphics(gcf,'Figs/sysID_poles.jpg',"Resolution",1000)


%% Gráfico en el tiempo
experimental = iddata(x_m(:,1),x_t(:,1),1/fs);
modelo = MFDID;
[y,fit,ic] = compare(experimental,modelo);
fit


gcf = figure('Position', [10 10 1000 300]);
subplot(2,1,1)
plot(y.SamplingInstants,x_m(:,1),'k','linewidth',1);
hold on;
plot(y.SamplingInstants,y.OutPutData,'r--','linewidth',1);
title('');
ylabel('Desp. [mm]');
xlabel('Tiempo [s]');
xlim([0 52])
legend('x_m experimental','x_m modelo','NumColumns',2,'Location','southeast')
grid on;

subplot(2,1,2)
plot(y.SamplingInstants,abs(x_m(:,1)-y.OutPutData),'k','linewidth',1)
legend(['NRMSE = ',num2str(goodnessOfFit(y.OutPutData,x_m(:,1),'NRMSE').*100),' %'])
xlabel('Tiempo [s]')
ylabel('|Error| [mm]')
xlim([0 52])
% ylim([0 1.5])
grid on

% exportgraphics(gcf,'Figs/sysID_T.jpg',"Resolution",1000)

%% Validación 
% Importar datos validación
addpath('Validacion');
Data = importdata('terremoto0x2A.txt').data;
Ts = 0.001; %s
x_t = (Data(20196:end,3)-Data(20196,3)).*10;
x_m = (Data(20196:end,1)-Data(20196,1)).*10;
experimental = iddata(x_m,x_t,Ts);
[~,~,~,delay] = Freq_Resp_Tong(x_t,x_m,1000); % Encontrar el retraso temporal

%% Figuras x_t vs x_m
modelo = MFDID;
[y,fit,ic] = compare(experimental,modelo);
x_m_model = y.OutPutData;
vector_tiempo = y.SamplingInstants;

gcf = figure('Position', [10 10 1000 300]);
subplot(2,4,[1,2,3])
plot(vector_tiempo,x_t,'k','linewidth',1);
hold on;
plot(vector_tiempo,x_m,'b--','linewidth',1);
title('');
ylabel('Desp. [mm]');
xlabel('Tiempo [s]');
xlim([0 65])
legend('x_t','x_m')
grid on;

subplot(2,4,[5,6,7])
plot(vector_tiempo,abs(x_t-x_m),'k','linewidth',1)
hold on;
plot(vector_tiempo,abs(x_t-x_m),'k','linewidth',1)
legend(['NRMSE = ',num2str(goodnessOfFit(x_m,x_t,'NRMSE').*100),' %'], ['Retraso = ',num2str(delay*1000),' ms'])
xlabel('Tiempo [s]')
ylabel('|Error| [mm]')
xlim([0 65])
grid on

subplot(2,4,[4,8])
plot(vector_tiempo,x_t,'k','linewidth',1);
hold on;
plot(vector_tiempo,x_m,'b--','linewidth',1);
title('');
ylabel('Desp. [mm]');
xlabel('Tiempo [s]');
xlim([45 47])
grid on;
exportgraphics(gcf,'Figs/sysID_xVal.jpg',"Resolution",1000)

%% x_m modelo vs x_m experimental
modelo = MFDID;
[y,fit,ic] = compare(experimental,modelo);
x_m_model = y.OutPutData;

gcf = figure('Position', [10 10 1000 300]);
subplot(2,4,[1,2,3])
plot(vector_tiempo,x_m,'k','linewidth',1);
hold on;
plot(vector_tiempo,x_m_model,'r--','linewidth',1);
title('');
ylabel('Desp. [mm]');
xlabel('Tiempo [s]');
xlim([0 65])
legend('x_m experimental','x_m modelo')
grid on;

subplot(2,4,[5,6,7])
plot(vector_tiempo,abs(x_m-x_m_model),'k','linewidth',1)
legend(['NRMSE = ',num2str(goodnessOfFit(x_m_model,x_m,'NRMSE').*100),' %'])
xlabel('Tiempo [s]')
ylabel('|Error| [mm]')
xlim([0 65])
% ylim([0 1.5])
grid on

subplot(2,4,[4,8])
plot(vector_tiempo,x_m,'k','linewidth',1);
hold on;
plot(vector_tiempo,x_m_model,'r--','linewidth',1);
title('');
ylabel('Desp. [mm]');
xlabel('Tiempo [s]');
xlim([45 47])
grid on;

% exportgraphics(gcf,'Figs/sysID_Val.jpg',"Resolution",1000)

