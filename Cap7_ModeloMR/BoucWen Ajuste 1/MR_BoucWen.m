function [fun, F_sim, F_lab, x, dx, t] = MR_BoucWen(Data,params)
% Entradas:
%     Data   : struct con los datos 
%     params : parametros del modelo

% Salidas:
%     fun    : valor del error RMSE
%     F_sim  : Fuerza obtenida con el modelo
%     F_lab  : Fuerza experimental
%     x      : desplazamiento experimental
%     dx     : velocidad experimental
%     t      : vector del tiempo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    options = simset('SrcWorkspace','current');
  
    t_data = Data.data(:,1);
    x_data = Data.data(:,2);
    dx_data = Data.data(:,3);
    F_data = Data.data(:,4);
    i = Data.Amp;
    Ts = t_data(2,1)-t_data(1,1);
    fs = 5000;           %Frecuencia de sampleo RTHS
    
    
    dtsim = 1/fs;
    totaltime = t_data(end,1);
    sim('AMR_ID_BoucWen.slx',[],options)
    t = F_AMR.Time;
    F_sim = F_AMR.Data;
    F_lab = F_lab.Data;
    x = x.Data;
    dx = dx.Data; 
     
    RMSE = sqrt(mean((F_lab - F_sim).^2))/max(F_lab);
    fun = RMSE;

    