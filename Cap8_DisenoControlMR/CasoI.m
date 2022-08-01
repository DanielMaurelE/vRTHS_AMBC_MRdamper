function [] = CasoI(terremoto,escala)

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

n = 2;   %order
fc = 100; %cutoff frequency
[numfilter,denfilter] = butter(n,fc/(fs/2)); %discrete filter coefficients

%% Modelo de la estructura
% Parámetros
m = 1000;  %masa por piso
d = 0.05;  %amortiguamiento por modo

Ms = eye(3)*m;  %Matriz masa
Ks = 10^7*[2.6055,-2.3134,0.5937;-2.3134,3.2561,-1.4420;0.5937,-1.4420,0.9267];  %Matriz rigidez benchmark

[mod,w2] = eig(Ks,Ms);
w = sqrt(w2);
f = w/2/pi*ones(3,1);

Cs = mod'^-1*(2*d*w)*mod^-1;        %Matriz amortiguamiento
n_pisos = length(Ms);
G1 = -Ms*ones(n_pisos,1);           %asociado al sismo
G2 = [-1;0;0];                      %asociado al control

% Modelo espacio-estado
ssA  = [zeros(n_pisos) eye(n_pisos);-Ms\Ks -Ms\Cs];
ssB1 = [zeros(n_pisos,1);Ms\G1];    % Entrada sismo
ssB2 = [zeros(n_pisos,1);Ms\G2];    % Entrada control
% salidas: [desp. rel., drifts, vel. rel., ace. abs]
ssC  = [eye(n_pisos) zeros(n_pisos);1 zeros(1,5);-1 1 zeros(1,4);0 -1 1 zeros(1,3);zeros(n_pisos) eye(n_pisos);-Ms\Ks -Ms\Cs];              
ssD1 = [zeros(4*n_pisos,1)];          % sismo
ssD2 = [zeros(3*n_pisos,1);Ms\G2];    % control 

Edificio = ss(ssA,[ssB1 ssB2],ssC,[ssD1 ssD2]);

%% Modelo Amortiguador MR
load Modelos/MR
params = parametros;

%% LQR
Edificio_LQR = ss(ssA,ssB2,ssC,ssD2);
R = 1;
Qdiag = [0,0,13,0,0,0];
Q = diag(10.^Qdiag);
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

%% Respuesta sistema 
Referencia = struct();
options = simset('SrcWorkspace','current');

tipo_control = 1; % 1: Passive-off, 2: Passive-on, 3: control LQG-COC
sim('Simulacion.slx',[],options)

% Tiempo
Referencia.t = y_sc.Time;

% Extraer datos sistema sin control
Referencia.sincontrol.y = y_sc.Data; 
Referencia.sincontrol.F = 0; 
Referencia.sincontrol.V = y_sc.Data(:,[10:12])*diag(Ms);

% Extraer datos en passive-off
Referencia.passiveoff.y = y_r.Data;
Referencia.passiveoff.F = F_r.Data;
Referencia.passiveoff.ic = ic_r.Data;
Referencia.passiveoff.V = y_r.Data(:,[10:12])*diag(Ms);

% Extraer datos en passive-on
tipo_control = 2; % 1: Passive-off, 2: Passive-on, 3: control LQG-COC
sim('Simulacion.slx',[],options)

Referencia.passiveon.y = y_r.Data;
Referencia.passiveon.F = F_r.Data;
Referencia.passiveon.ic = ic_r.Data;
Referencia.passiveon.V = y_r.Data(:,[10:12])*diag(Ms);


% Extraer datos en LQG-COC
tipo_control = 3; % 1: Passive-off, 2: Passive-on, 3: control LQG-COC
sim('Simulacion.slx',[],options)

Referencia.LQGCOC.y = y_r.Data;
Referencia.LQGCOC.x_kalman = x_kalman.Data;
Referencia.LQGCOC.F = F_r.Data;
Referencia.LQGCOC.ic = ic_r.Data;
Referencia.LQGCOC.V = y_r.Data(:,[10:12])*diag(Ms);


%% Sistema Activo ideal
sim('Activo_ideal.slx',[],options)
Referencia.Activo.t = y_A.Time;
Referencia.Activo.y = y_A.Data;
Referencia.Activo.F = F_A.Data;
Referencia.Activo.V = y_A.Data(:,[10:12])*diag(Ms);


%% Guardar resultados
save(['CasoI_Resultados/',terremoto,'/',terremoto,num2str(escala)],'Referencia')