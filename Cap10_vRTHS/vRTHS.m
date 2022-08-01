function [] = vRTHS(terremoto,escala,caso)
options = simset('SrcWorkspace','current');

%% Terremoto de entrada
switch terremoto
    case 'Centro'
        Ug = load('Terremotos/ElCentroAccelNoScaling.mat').ElCentroAccel';       
    case 'Kobe'
        Ug = load('Terremotos/KobeAccelNoScaling.mat').KobeAccel';
        
    case 'Maule' 
        Ug = load('Terremotos/MauleAccelNoScaling.mat').MauleAccel';
    otherwise 
        error('Terremoto desconocido')
end

Ug(:,2) = Ug(:,2)*escala/100;

%% Parámetros de la simulación
% Tiempo de muestreo y total:
fs = 2^12;          % Frecuencia de muestreo  
dtsim = 1/fs;       % Tiempo de muestreo 
ttotal = Ug(end,1); % Tiempo total de la simulación

% Sensores:
% Sensor de desplazamiento
rmsdesp = 2e-5;  % rms m
rms_noise_desp = rmsdesp^2*dtsim;    %  Noise power

% Sensor de fuerza
rmsF = 20;  %rms N
rms_noise_F = rmsF^2*dtsim; % Noise power

% Límites de saturación 
sat_limit_upper = +3.8;         % Volts
sat_limit_lower = -3.8;         % Volts

% Intervalo de cuantización
quantize_int = 1 / 2^18;        % 18 bit channel

%% Modelo de la estructura
% Parámetros
switch caso
    case 'CasoI'
        m = 1000;  %masa por piso
        d = 0.05;  %amortiguamiento por modo
    case 'CasoII'
        m = 1100;  %masa por piso
        d = 0.04;  %amortiguamiento por modo    
    case 'CasoIII' 
        m = 1300;  %masa por piso
        d = 0.03;  %amortiguamiento por modo
    otherwise 
        error('Caso desconocido')
end

Ms = eye(3)*m;  %Matriz masa
Ks = 10^7*[2.6055,-2.3134,0.5937;-2.3134,3.2561,-1.4420;0.5937,-1.4420,0.9267];  %Matriz rigidez benchmark

[mod,w2] = eig(Ks,Ms);
w = sqrt(w2);
f = w/2/pi*ones(3,1);

Cs = mod'^-1*(2*d*w)*mod^-1;        %Matriz amortiguamiento
n_pisos = length(Ms);
G1 = -Ms*ones(n_pisos,1);           %asociado al sismo
G2 = zeros(n_pisos,1); G2(1) = -1;  %asociado al control

% Modelo espacio-estado
ssA  = [zeros(n_pisos) eye(n_pisos);-Ms\Ks -Ms\Cs];
ssB1 = [zeros(n_pisos,1);Ms\G1];    % Entrada sismo
ssB2 = [zeros(n_pisos,1);Ms\G2];    % Entrada control
ssC  = [eye(2*n_pisos);-Ms\Ks -Ms\Cs];             % Aceleraciones absolutas 
ssD1 = [zeros(2*n_pisos,1);zeros(n_pisos,1)];          % sismo
ssD2 = [zeros(2*n_pisos,1);Ms\G2];                     % control 

Edificio = ss(ssA,[ssB1 ssB2],ssC,[ssD1 ssD2]);

%% Modelo Actuador
load Modelos/ACTUADOR
G0 = idtf(sys_mfdid);
[act_num,act_den] = tfdata(G0,'v');
act_num = act_num(end);

%% Compensador AMB
A_amb_i = flip(act_den/act_num);        %Initial compensator's parameters a_i
A_amb_i = [A_amb_i 0];

%Noise filter for adaptation process
n = 4;   %order
fc = 20; %cutoff frequency
[numfilter,denfilter] = butter(n,fc/(fs/2)); %discrete filter coefficients

% Adaptive gains
op_location = load('Modelos/OPTIMO').op_location';
adaptivegain = diag(10.^op_location);   %obtained from calibration

% Backward difference
u=[1 0 0 0];
dudt=[1 -1 0 0]/dtsim;
dudt2=[1 -2 1 0]/dtsim^2;
dudt3=[1 -3 3 -1]/dtsim^3;
Derivadas=[u;dudt;dudt2;dudt3];  %FIR filter for numerical derivatives

%Adaptive parameters constraints
%a0
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

%% Modelo Amortiguador MR
load Modelos/MR
params = parametros;

%% LQR
Edificio_LQR = ss(ssA,ssB2,ssC,ssD2);
Q = diag([0,0,1e13,0,0,0]);
R = 1;
[Kc,S,E] = lqr(Edificio_LQR,Q,R);

%% Filtro de Kalman
% Covarianzas
Q_w = 0.1;
R_v = 5e-4;

% Modelo espacio-estado
ssA  = [zeros(n_pisos) eye(n_pisos);-Ms\Ks -Ms\Cs];
ssB1 = [zeros(n_pisos,1);Ms\G1];    % Entrada sismo
ssB2 = [zeros(n_pisos,1);Ms\G2];    % Entrada control
ssC  = [-Ms\Ks -Ms\Cs];             % Aceleraciones absolutas 
ssD1 = [zeros(n_pisos,1)];          % sismo
ssD2 = [Ms\G2];                     % control 

Edificio_kalman = ss(ssA,[ssB2 ssB1],ssC,[ssD2 ssD1]);
sys_kalman = kalman(Edificio_kalman,Q_w,R_v);
sys_kalman = sys_kalman(4:9,:);

%% Respuesta sistema de referencia
Referencia = struct();
SinControl = struct();

tipo_control = 1; % 1: Passive-off, 2: Passive-on, 3: control LQG-COC
sim('vRTHS_Referencia.slx',[],options)

% Tiempo
SinControl.t = y_sc.Time;
Referencia.t = y_r.Time;

% Extraer datos sistema sin control
SinControl.y = y_sc.Data; 

% Extraer datos en passive-off
Referencia.passiveoff.y = y_r.Data;
Referencia.passiveoff.F = F_r.Data;
Referencia.passiveoff.ic = ic_r.Data;

% Extraer datos en passive-on
tipo_control = 2; % 1: Passive-off, 2: Passive-on, 3: control LQG-COC
sim('vRTHS_Referencia.slx',[],options)

Referencia.passiveon.y = y_r.Data;
Referencia.passiveon.F = F_r.Data;
Referencia.passiveon.ic = ic_r.Data;

% Extraer datos en LQG-COC
tipo_control = 3; % 1: Passive-off, 2: Passive-on, 3: control LQG-COC
sim('vRTHS_Referencia.slx',[],options)

Referencia.LQGCOC.y = y_r.Data;
Referencia.LQGCOC.x_kalman = x_kalman.Data;
Referencia.LQGCOC.F = F_r.Data;
Referencia.LQGCOC.ic = ic_r.Data;

%% Guardar resultados
save(['Resultados/',caso,'/',terremoto,num2str(escala),'SinControl'],'SinControl')
save(['Resultados/',caso,'/',terremoto,num2str(escala),'Referencia'],'Referencia')

%% Respuesta sistemas vRTHS
SinCompensar = struct();
ConCompensar = struct();

%% Passive-off
tipo_control = 1; % 1: Passive-off, 2: Passive-on, 3: control LQG-COC
sim('vRTHS_Compensacion.slx',[],options)

% Tiempo
SinCompensar.t = y0.Time;

ConCompensar.t = y1.Time;

% Desplazamientos, velocidades y aceleraciones Sub. numérica 
SinCompensar.passiveoff.y = y0.Data;

ConCompensar.passiveoff.y = y1.Data;

% Desplazamiento medido(m), comandado (c) y objetivo (t)
SinCompensar.passiveoff.x_m  = x_m0.Data;
SinCompensar.passiveoff.x_t  = x_t0.Data;

ConCompensar.passiveoff.x_m  = x_m1.Data;
ConCompensar.passiveoff.x_t  = x_t1.Data;
ConCompensar.passiveoff.x_c  = x_c1.Data;

% Resultados de fuerza
SinCompensar.passiveoff.F  = F0.Data;  

ConCompensar.passiveoff.F  = F1.Data;  

% Resultado de corriente aplicada
SinCompensar.passiveoff.ic = ic0.Data;

ConCompensar.passiveoff.ic = ic1.Data;

% Parámetros adaptivos
ConCompensar.passiveoff.amb = amb.Data;

%% Passive-on
tipo_control = 2; % 1: Passive-off, 2: Passive-on, 3: control LQG-COC
sim('vRTHS_Compensacion.slx',[],options)

% Desplazamientos, velocidades y aceleraciones Sub. numérica 
SinCompensar.passiveon.y = y0.Data;

ConCompensar.passiveon.y = y1.Data;

% Desplazamiento medido(m), comandado (c) y objetivo (t)
SinCompensar.passiveon.x_m  = x_m0.Data;
SinCompensar.passiveon.x_t  = x_t0.Data;

ConCompensar.passiveon.x_m  = x_m1.Data;
ConCompensar.passiveon.x_t  = x_t1.Data;
ConCompensar.passiveon.x_c  = x_c1.Data;

% Resultados de fuerza
SinCompensar.passiveon.F  = F0.Data;  

ConCompensar.passiveon.F  = F1.Data;  

% Resultado de corriente aplicada
SinCompensar.passiveon.ic = ic0.Data;

ConCompensar.passiveon.ic = ic1.Data;

% Parámetros adaptivos
ConCompensar.passiveon.amb = amb.Data;

%% LQG-COC
tipo_control = 3; % 1: Passive-off, 2: Passive-on, 3: control LQG-COC
sim('vRTHS_Compensacion.slx',[],options)

% Desplazamientos, velocidades y aceleraciones Sub. numérica 
SinCompensar.LQGCOC.y = y0.Data;

ConCompensar.LQGCOC.y = y1.Data;

% Desplazamientos, velocidades y aceleraciones estimados por F. Kalman
SinCompensar.LQGCOC.x_kalman = x_kalman0.Data;

ConCompensar.LQGCOC.x_kalman = x_kalman1.Data;

% Desplazamiento medido(m), comandado (c) y objetivo (t)
SinCompensar.LQGCOC.x_m  = x_m0.Data;
SinCompensar.LQGCOC.x_t  = x_t0.Data;

ConCompensar.LQGCOC.x_m  = x_m1.Data;
ConCompensar.LQGCOC.x_t  = x_t1.Data;
ConCompensar.LQGCOC.x_c  = x_c1.Data;

% Resultados de fuerza
SinCompensar.LQGCOC.F  = F0.Data;  

ConCompensar.LQGCOC.F  = F1.Data;  

% Resultado de corriente aplicada
SinCompensar.LQGCOC.ic = ic0.Data;

ConCompensar.LQGCOC.ic = ic1.Data;

% Parámetros adaptivos
ConCompensar.LQGCOC.amb = amb.Data;

%% Guardar resultados
save(['Resultados/',caso,'/',terremoto,num2str(escala),'SinCompensar'],'SinCompensar')
save(['Resultados/',caso,'/',terremoto,num2str(escala),'ConCompensar'],'ConCompensar')

