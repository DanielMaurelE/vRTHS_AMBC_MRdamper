function in = indicadores(SinControl,Referencia,Datos,graf,compensado)
if nargin < 4
    graf = 0;
end
if nargin < 5
    compensado = 0;
end

t = SinControl.t;
x_m = Datos.x_m;
x_t = Datos.x_t;
x_r = Referencia.y(:,1);
F_r = Referencia.F;
F_m = Datos.F;
y3_sc = SinControl.y(:,3);
y3 = Datos.y(:,3);
ddy3_sc = SinControl.y(:,9);
ddy3 = Datos.y(:,9);

J1 = goodnessOfFit(x_m,x_r,'NRMSE')*100;
J2 = goodnessOfFit(F_m,F_r,'NRMSE')*100;
J3 = goodnessOfFit(x_m,x_t,'NRMSE')*100;
J4 = (max(abs(y3)))./max(abs(y3_sc))*100;
J5 = (max(abs(ddy3)))./max(abs(ddy3_sc))*100;
J6 = max(abs(F_m))./1000;
J7 = (rms(y3))./rms(y3_sc)*100;
J8 = (rms(ddy3))./rms(ddy3_sc)*100;
J9 = rms(F_m)./1000;
in = [J1,J2,J3,J4,J5,J6,J7,J8,J9];

%% Gráficos
if compensado == 1
        aux = 'compensado';
    else
        aux = 'no_compensado';
end
    
if graf == 1
    % J1 y J2
    % x_m vs x_r
    gcf = figure('Position', [10 10 1000 300]);
    tiledlayout(2,1);
    hf = nexttile; 
    plot(t,x_r*1000,'k','linewidth',1,'DisplayName','Sistema de referencia');
    hold on;
    plot(t,x_m*1000,'r--','linewidth',1,'DisplayName','vRTHS');
    ylim([-0.6 0.6])
    title('');
    ylabel('$x_1$ [mm]','Interpreter','Latex','FontSize',14)
    xlabel('Tiempo [s]','Interpreter','Latex');
    xlim([20 30])
    grid on;
    title(['J_1 = ',num2str(J1,3),' %'])
    lh = legend('Orientation','horizontal');
%     lh.Layout.Tile = 'South';

    % F_m vs F_r
    hf = nexttile; 
    plot(t,F_r./1000,'k','linewidth',1,'DisplayName','Sistema de referencia');
    hold on;
    plot(t,F_m./1000,'r--','linewidth',1,'DisplayName','vRTHS');
    title('');
    ylim([-5 5])
    ylabel('F [kN]','Interpreter','Latex','FontSize',14)
    xlabel('Tiempo [s]')
    xlim([20 30])
    grid on;
    ax = gca;
    title(['J_2 = ',num2str(J2,3),' %'],'Interpreter','Latex')


    exportgraphics(gcf,strcat('Figs/J1-J2_',aux,'.jpg'),"Resolution",1000)
    
    % J3 
    % x_m vs x_t
    gcf = figure('Position', [10 10 1000 150]);
    tiledlayout(1,1);
    hf = nexttile; 
    plot(t,x_t*1000,'b','linewidth',1,'DisplayName','x_t');
    hold on;
    plot(t,x_m*1000,'r--','linewidth',1,'DisplayName','x_m');
    title('');
    ylim([-1 1])
    ylabel('$x_1$ [mm]','Interpreter','Latex','FontSize',14)
    xlabel('Tiempo [s]');
    xlim([20 30])
    grid on;
    title(['J_3 = ',num2str(J3,3),' %'])
    lh = legend('Orientation','horizontal');
%     lh.Layout.Tile = 'South';
    
    exportgraphics(gcf,strcat('Figs/J3',aux,'.jpg'),"Resolution",1000)

    % J4, J5, J6, J7, J8, J9 
    % desplazamiento
    gcf = figure('Position', [10 10 1000 500]);
    tiledlayout(3,1);
    hf = nexttile; 
    plot(t,y3_sc.*1000,'k','linewidth',1,'DisplayName','Sistema sin control');
    hold on;
    plot(t,y3.*1000,'r--','linewidth',1,'DisplayName','vRTHS');
    title('');
    ylim([-10 10])
    ylabel('$x_3$ [mm]','Interpreter','Latex','FontSize',14)
    xlabel('Tiempo [s]');
    xlim([5 35])
    grid on;
    title(['J_{4} = ',num2str(J4,3),' %, J_{7} = ',num2str(J7,3),' %'])
    lh = legend('Orientation','horizontal');

    % aceleracion 
    hf = nexttile; 
    plot(t,ddy3_sc.*1000,'k','linewidth',1,'DisplayName','Sistema sin control');
    hold on;
    plot(t,ddy3.*1000,'r--','linewidth',1,'DisplayName','vRTHS');
    title('');
    ylim([-10000 10000])
    ylabel('$\ddot{x}_3$ [mm/s$^2$]','Interpreter','Latex','FontSize',14)
    xlabel('Tiempo [s]');
    xlim([5 35])
    grid on;
    ax = gca;
    title(['J_{5} = ',num2str(J5,3),' %, J_{8} = ',num2str(J8,3),' %'])
    
    
    % F_m vs F_r
    hf = nexttile; 
    plot(t,F_m./1000,'r','linewidth',1,'DisplayName','vRTHS');
    title('');
    ylim([-5 5])
    ylabel('F [kN]','Interpreter','Latex','FontSize',14)
    xlabel('Tiempo [s]');
    xlim([5 35])
    grid on;
    ax = gca;
    title(['J_{6} = ',num2str(J6,3),' kN, J_{9} = ',num2str(J9,3),' kN'])
    
    
    exportgraphics(gcf,strcat('Figs/J4-J9_',aux,'.jpg'),"Resolution",1000)


    %% Parámetros AMB 
    if isfield(Datos,'amb')
        amb = Datos.amb;
        totaltime = t(end);

        gcf = figure('Position', [10 10 500 300]);
        tiledlayout(2,2);
        hf = nexttile; 
        plot(t,amb(:,1),'k')
        xlabel('Tiempo [s]')
        ylabel('a_0')
        xlim([0 totaltime])
        grid on

        hf = nexttile; 
        plot(t,amb(:,2),'k')
        xlabel('Tiempo [s]')
        ylabel('a_1  [s]')
        xlim([0 totaltime])
        grid on

        hf = nexttile; 
        plot(t,amb(:,3),'k')
        xlabel('Tiempo [s]')
        ylabel('a_2  [s^2]')
        xlim([0 totaltime])
        grid on

        hf = nexttile; 
        plot(t,amb(:,4),'k')
        xlabel('Tiempo [s]')
        ylabel('a_3   [s^3]')
        xlim([0 totaltime])
        grid on
        
        exportgraphics(gcf,strcat('Figs/adaptivegains.jpg'),"Resolution",1000)
    end
end


