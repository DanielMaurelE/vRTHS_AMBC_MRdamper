%% Adaptive gains calibration 

%% Calibration plants
%Actuator parameters
%% Transfer system dynamics and Initial Model
load modelo
G0 = idtf(sys_mfdid);
[mag,phase,wout0] = bode(G0,[0.1:0.1:15*2*pi]);
mag = squeeze(mag);
phase = squeeze(phase);
wout0 = wout0./2/pi;

[num,den]=tfdata(G0,'v');
num = num(end);
GpLoadSysred=tf(num,den);    %Initial model
AMB=tf(den/num,1);           
A_amb_i = [flip(den/num) 0]        %Initial compensator's parameters a_i

nsimu=50;   %Number of realization for visualization of random control plants
fbode=0.1:0.1:15;  %bode frecuency range
wbode=2*pi*fbode;
[maginicial,phaseinicial,wout] = bode(tf(1,flip(A_amb_i)),wbode);

b1cals = [0.96, 1.01];
b2cals = [0.02, 0.03];
b3cals = [1e-5, 3.5e-4];

%% Calibration earthquake
load('ElCentroAccelNoScaling.mat');

Ugcal=ElCentroAccel';

dtUg = Ugcal(2,1)-Ugcal(1,1);   %record time step
xUg = detrend(cumtrapz(cumtrapz(Ugcal(:,2)))*dtUg^2).*1000; % record displacement, mm

%Portion of the earthquake utilized for calibration to improve efficiency
Ugdesde = 5;         %from sec
Ughasta = 12;        %to sec
Ugcal=[(0:dtUg:(Ughasta-Ugdesde))',Ugcal(Ugdesde/dtUg+1:Ughasta/dtUg+1,2)];

tcal=max(Ugcal(:,1));  %total time of calibration simulation

escalas=[30,60]/100;  %Scaling factor range for calibration (uniform distribution)

%% SDOF Numerical substructures for calibration

frefcals=[2.8,4];   %natural frequency range (Hz)
drefcals=[3,5]/100; %damping ratio range

%% Matlab function for calibration

%input: adaptive gain matrix 
%output: mean J2 error (R2) , std of J2 (dR2) , max J2 (maxR2)

%number of realizations
realizations = 100;


%function (x=log10(diag(adaptive gain matrix))
fun=@(x) AMB_6_R2function(x,escalas,frefcals,drefcals,realizations,b1cals,b2cals,b3cals);

%example with Gamma=diag(10.^[7,4,1,-2])
tic
% [prom,desv,maxj2] = fun([7.3,4.8,2.7,-2.5])
toc



