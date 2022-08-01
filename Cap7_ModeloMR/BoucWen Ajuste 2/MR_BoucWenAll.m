function fun = MR_BoucWenAll(Datas,param)
% Entradas:
%     Data   : struct con los datos 
%     param  : parametros del modelo

% Salidas:
%     fun    : valor de la suma de los RMSE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    options = simset('SrcWorkspace','current');
    fun = 0;
    
    for Data = Datas
        t_data = Data.data(:,1);
        x_data = Data.data(:,2);
        dx_data = Data.data(:,3);
        F_data = Data.data(:,4);
        i = Data.Amp;
        Ts = t_data(2,1)-t_data(1,1);
        fs = 5000;           %Frecuencia de sampleo RTHS 

        k0 = param(1);
        k1 = param(2);
        x0 = param(3);
        alpha = param(4)*i^2+param(5)*i+param(6);
        gamma = param(7)*exp(param(8)*i);
        c0 = param(9)*i+param(10);
        c1 = param(11)*i+param(12);
        params = [k0,k1,x0,alpha,gamma,c0,c1];
        
        dtsim = 1/fs;
        totaltime = t_data(end,1);
        sim('AMR_ID_BoucWen.slx',[],options)
        t = F_AMR.Time;
        F_sim = F_AMR.Data;
        F_lab = F_lab.Data;
        x = x.Data;
        dx = dx.Data; 

        zci = @(v) find(diff(sign(v))); % funcion para encontrar los cruces por 0

        pos = zci(F_lab);
        index = floor(length(pos)/2);
        F_sim = F_sim(pos(index):pos(index+2));
        F_lab = F_lab(pos(index):pos(index+2));
        x = x(pos(index):pos(index+2));
        dx = dx(pos(index):pos(index+2));
        t = t(pos(index):pos(index+2));
        t = t-t(1);

        RMSE = sqrt(mean((F_lab - F_sim).^2))/max(F_lab);
        fun = fun+RMSE;
    end

    