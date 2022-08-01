%% Definición parámetros simulación
% María Quiroz
% Estudiante de Magíster en Ciencias de la Ingeniería Civil, UTFSM
clear all; close all; clc

%% Definicion 
%Simulation parameters
fs=4096;           %Sampling frequency   
dtsim=1/fs;        %Time step

% Sensores:
%Displacement sensor
rms_noise_desp = 6.25e-14;    %  Noise power
rmsdesp = sqrt(rms_noise_desp/dtsim)*1000;  %rms mm

%Force sensor
rms_noise_F = 1.16e-3; % Noise power
rmsF = sqrt(rms_noise_F/dtsim);  %rms N

% Saturation limits 
sat_limit_upper = +3.8;         % Volts
sat_limit_lower = -3.8;         % Volts
% Quantization interval
quantize_int = 1 / 2^18;        % 18 bit channel

%% Transfer system dynamics and Initial Model
load modelo
G0 = idtf(sys_mfdid);
[mag,phase,wout0] = bode(G0,[0.1:0.1:20*2*pi]);
mag = squeeze(mag);
phase = squeeze(phase);
wout0 = wout0./2/pi;

[num,den]=tfdata(G0,'v');
num = num(end);
GpLoadSysred=tf(num,den);    %Initial model
AMB=tf(den/num,1);           
A_amb_i = flip(den/num)        %Initial compensator's parameters a_i

%% Compensator

% Backward difference
u=[1 0 0 0];
dudt=[1 -1 0 0]/dtsim;
dudt2=[1 -2 1 0]/dtsim^2;
dudt3=[1 -3 3 -1]/dtsim^3;
Derivadas=[u;dudt;dudt2;dudt3];  %FIR filter for numerical derivatives

%Adaptive parameters constraints
%a0
amb0max=inf;   
amb0min=0;
%a1
amb1max=inf;
amb1min=0;
%a2
amb2max=inf;
amb2min=0;
%a3
amb3max=inf;
amb3min=0;
%For saturation block
max_amb=[amb0max,amb1max,amb2max,amb3max];
min_amb=[amb0min,amb1min,amb2min,amb3min];

%Noise filter for adaptation process
n = 4;   %order
fc = 20; %cutoff frequency
% fs = 1000;
[numfilter,denfilter] = butter(n,fc/(fs/2)); %discrete filter coefficients

G_filter = tf(numfilter,denfilter,1/fs);
[mag_filter,phase_filter,w] = bode(G_filter,[0.1:0.1:20*2*pi]);
mag_filter = squeeze(mag_filter);
phase_filter = squeeze(phase_filter);
w = w./2/pi;

mags = []; phases = [];
b1cals = [0.96, 1.01];
b2cals = [0.02, 0.03];
b3cals = [1e-5, 3.5e-4];
latin = lhsdesign(100,3); %latin hypercube sampling
b1par = b1cals(1)+(b1cals(2)-b1cals(1))*latin(:,1);
b2par = b2cals(1)+(b2cals(2)-b2cals(1))*latin(:,2);
b3par = b3cals(1)+(b3cals(2)-b3cals(1))*latin(:,3);

for i = 1:100
    den = [b3par(i),b2par(i),b1par(i)];
    G = tf(1,den);
    [m,p,~] = bode(G,[0.1:0.1:20*2*pi]);
    mags(:,i) = squeeze(m);
    phases(:,i) = wrapTo180(squeeze(p));
end
mags_max = db(max(mags,[],2));
mags_min = db(min(mags,[],2));
phases_max = max(phases,[],2);
phases_min = min(phases,[],2);
delays = wrapTo180(phases)./360./w.*1000;
delays_max = max(delays,[],2);
delays_min = min(delays,[],2);

gcf = figure('Position', [10 10 800 400]);
subplot(3,1,1);
plot(wout0,db(abs(mag)),'r','Linewidth',2); hold on;
plot(w,db(abs(mag_filter)),'g--','Linewidth',2);
patch([w' fliplr(w')], [mags_min' fliplr(mags_max')], 'k' ,'FaceAlpha',.2,'EdgeColor','none')
% plot(w,db(abs(mags)),'Color',grayColor,'Linewidth',1);
xlim([0 20])
grid on;
legend('Modelo inicial','Filtro Butterworth','Plantas de calibración','Location','southwest')
ylabel('Magnitud [dB]');

subplot(3,1,2);
plot(wout0,wrapTo180(phase),'r','Linewidth',2'); hold on;
plot(w,wrapTo180(phase_filter),'g--','Linewidth',2');
patch([w' fliplr(w')], [phases_min' fliplr(phases_max')], 'k' ,'FaceAlpha',.2,'EdgeColor','none')
xlim([0 20])
grid on;
ylabel('Fase [°]');

subplot(3,1,3);
plot(wout0,phase./360./wout0.*1000,'r','Linewidth',2); hold on;
plot(w,wrapTo180(phase_filter)./360./w.*1000,'g--','Linewidth',2);
patch([w' fliplr(w')], [delays_min' fliplr(delays_max')], 'k' ,'FaceAlpha',.2,'EdgeColor','none')
xlim([0 20])
% ylim([-200 0])
grid on;
ylabel('Retraso [ms]');
xlabel('Frecuencia [Hz]');
% exportgraphics(gcf,'Figs/act_initialmodel.jpg',"Resolution",1000)


