
%% Cargar datos
Centro = load(['CasoI_Resultados/Centro/Centro50'],'Referencia').Referencia;
Kobe = load(['CasoI_Resultados/Kobe/Kobe50'],'Referencia').Referencia;
Maule = load(['CasoI_Resultados/Maule/Maule50'],'Referencia').Referencia;

%% Corte basal
gcf = figure('Position', [10 10 1000 500]);
tiledlayout(3,1);

Vsc = Centro.sincontrol.V;
Vcc = Centro.LQGCOC.V;
Vmax = (max(abs(Vsc))-max(abs(Vcc)))./max(abs(Vsc))*100;
t = Centro.t;
hf = nexttile; 
plot(t,Vsc./1000,'k','linewidth',1,'DisplayName','Sin control'); 
hold on;
plot(t,Vcc./1000,'r--','linewidth',1,'DisplayName','Semi-activo'); 
xlabel('Tiempo [s]');
ylabel('Corte basal [kN]')
grid on;
xlim([0 t(end)]);
title('Terremoto El Centro 1940');
text(hf.XLim(end)*0.6,hf.YLim(end)*0.85,['Reducción del maximo de ',num2str(Vmax,2),' %'],'FontSize',8)

Vsc = Kobe.sincontrol.V;
Vcc = Kobe.LQGCOC.V;
t = Kobe.t;
Vmax = (max(abs(Vsc))-max(abs(Vcc)))./max(abs(Vsc))*100;
hf = nexttile; 
plot(t,Vsc./1000,'k','linewidth',1,'DisplayName','Sin control'); 
hold on;
plot(t,Vcc./1000,'r--','linewidth',1,'DisplayName','Semi-activo'); 
xlabel('Tiempo [s]');
ylabel('Corte basal [kN]')
grid on;
xlim([0 t(end)]);
title('Terremoto Kobe 1995');
text(hf.XLim(end)*0.6,hf.YLim(end)*0.85,['Reducción del maximo de ',num2str(Vmax,2),' %'],'FontSize',8)


Vsc = Maule.sincontrol.V;
Vcc = Maule.LQGCOC.V;
t = Maule.t;
Vmax = (max(abs(Vsc))-max(abs(Vcc)))./max(abs(Vsc))*100;

hf = nexttile; 
plot(t,Vsc./1000,'k','linewidth',1,'DisplayName','Sin control'); 
hold on;
plot(t,Vcc./1000,'r--','linewidth',1,'DisplayName','Semi-activo'); 
xlabel('Tiempo [s]');
ylabel('Corte basal [kN]')
grid on;
xlim([0 t(end)]);
text(hf.XLim(end)*0.6,hf.YLim(end)*0.85,['Reducción del maximo de ',num2str(Vmax,2),' %'],'FontSize',8)
title('Terremoto Maule 2010');
lh = legend('Orientation','horizontal');
lh.Layout.Tile = 'South';
% exportgraphics(gcf,'Figs/CasoI/Corte_casoI.jpg',"Resolution",1000)

%% Gráficos max
index = 1;
gcf = figure('Position', [10 10 600 800]);
tiledlayout(3,3);
terremotos = {'El Centro 50%','Kobe 50%','Maule 50%'};
for i = [Centro, Kobe, Maule]
    Referencia = i;
    % Desplazamiento máx por piso
    pisos = [1,2,3];

    nexttile 
    plot(max(abs(Referencia.sincontrol.y(:,1:3)*1000)),pisos,'-ok','linewidth',2,'MarkerFaceColor','k'); 
    hold on;
    plot(max(abs(Referencia.passiveoff.y(:,1:3)*1000)),pisos,'-x','linewidth',1.5); 
    plot(max(abs(Referencia.passiveon.y(:,1:3)*1000)),pisos,'-*','linewidth',1.5); 
    plot(max(abs(Referencia.LQGCOC.y(:,1:3)*1000)),pisos,'--o','linewidth',1.5); 
    plot(max(abs(Referencia.Activo.y(:,1:3)*1000)),pisos,'--p','linewidth',1.5); 
    xlabel('$x_{RMS}$ [mm]','Interpreter','Latex')
    ylabel('Piso')
    grid on;
    xlim([0 20])
    yticks(pisos)
    xlm=xlim;
    title(terremotos(index),'Position',[xlm(2)/2 xlm(1)+3.05]);

    nexttile 
    plot(max(abs(Referencia.sincontrol.y(:,4:6)*1000)),pisos,'-ok','linewidth',2,'MarkerFaceColor','k'); 
    hold on;
    plot(max(abs(Referencia.passiveoff.y(:,4:6)*1000)),pisos,'-x','linewidth',1.5); 
    plot(max(abs(Referencia.passiveon.y(:,4:6)*1000)),pisos,'-*','linewidth',1.5); 
    plot(max(abs(Referencia.LQGCOC.y(:,4:6)*1000)),pisos,'--o','linewidth',1.5);
    plot(max(abs(Referencia.Activo.y(:,4:6)*1000)),pisos,'--p','linewidth',1.5);
    xlabel('$d_{RMS}$ [mm]','Interpreter','Latex')
    ylabel('Piso')
    grid on;
    xlim([0 10])
    yticks(pisos)
    xlm=xlim;
    title(terremotos(index),'Position',[xlm(2)/2 xlm(1)+3.05]);

    nexttile 
    plot(max(abs(Referencia.sincontrol.y(:,10:12))),pisos,'-ok','DisplayName','Sin control','linewidth',2,'MarkerFaceColor','k'); 
    hold on;
    plot(max(abs(Referencia.passiveoff.y(:,10:12))),pisos,'-x','DisplayName','Passive-off','linewidth',1.5); 
    plot(max(abs(Referencia.passiveon.y(:,10:12))),pisos,'-*','DisplayName','Passive-on','linewidth',1.5); 
    plot(max(abs(Referencia.LQGCOC.y(:,10:12))),pisos,'--o','DisplayName','LQG-COC','linewidth',1.5); 
    plot(max(abs(Referencia.Activo.y(:,10:12))),pisos,'--p','DisplayName','Activo','linewidth',1.5);
    xlabel('$\ddot{x}_{RMS}$ [m/s$^2$]','Interpreter','Latex')
    ylabel('Piso')
    grid on;
    xlim([0 9])
    yticks(pisos)
    xlm=xlim;
    title(terremotos(index),'Position',[xlm(2)/2 xlm(1)+3.05]);
    index = index+1;
end
lh = legend('Orientation','horizontal');
lh.Layout.Tile = 'South';
% exportgraphics(gcf,['Figs/CasoI/RespuestaMAX_casoI.jpg'],"Resolution",1000)

%% Gráficos RMS
index = 1;
gcf = figure('Position', [10 10 600 800]);
tiledlayout(3,3);
terremotos = {'El Centro 50%','Kobe 50%','Maule 50%'};
for i = [Centro, Kobe, Maule]
    Referencia = i;
    % Desplazamiento máx por piso
    pisos = [1,2,3];

    nexttile 
    plot(rms(Referencia.sincontrol.y(:,1:3)*1000),pisos,'-ok','linewidth',2,'MarkerFaceColor','k'); 
    hold on;
    plot(rms(Referencia.passiveoff.y(:,1:3)*1000),pisos,'-x','linewidth',1.5); 
    plot(rms(Referencia.passiveon.y(:,1:3)*1000),pisos,'-*','linewidth',1.5); 
    plot(rms(Referencia.LQGCOC.y(:,1:3)*1000),pisos,'--o','linewidth',1.5); 
    plot(rms(Referencia.Activo.y(:,1:3)*1000),pisos,'--p','linewidth',1.5); 
    xlabel('$x_{RMS}$ [mm]','Interpreter','Latex')
    ylabel('Piso')
    grid on;
    xlim([0 3])
    yticks(pisos)
    xlm=xlim;
    title(terremotos(index),'Position',[xlm(2)/2 xlm(1)+3.05]);

    nexttile 
    plot(rms(Referencia.sincontrol.y(:,4:6)*1000),pisos,'-ok','linewidth',2,'MarkerFaceColor','k'); 
    hold on;
    plot(rms(Referencia.passiveoff.y(:,4:6)*1000),pisos,'-x','linewidth',1.5); 
    plot(rms(Referencia.passiveon.y(:,4:6)*1000),pisos,'-*','linewidth',1.5); 
    plot(rms(Referencia.LQGCOC.y(:,4:6)*1000),pisos,'--o','linewidth',1.5);
    plot(rms(Referencia.Activo.y(:,4:6)*1000),pisos,'--p','linewidth',1.5);
    xlabel('$d_{RMS}$ [mm]','Interpreter','Latex')
    ylabel('Piso')
    grid on;
    xlim([0 1.5])
    yticks(pisos)
    xlm=xlim;
    title(terremotos(index),'Position',[xlm(2)/2 xlm(1)+3.05]);

    nexttile 
    plot(rms(Referencia.sincontrol.y(:,10:12)),pisos,'-ok','DisplayName','Sin control','linewidth',2,'MarkerFaceColor','k'); 
    hold on;
    plot(rms(Referencia.passiveoff.y(:,10:12)),pisos,'-x','DisplayName','Passive-off','linewidth',1.5); 
    plot(rms(Referencia.passiveon.y(:,10:12)),pisos,'-*','DisplayName','Passive-on','linewidth',1.5); 
    plot(rms(Referencia.LQGCOC.y(:,10:12)),pisos,'--o','DisplayName','LQG-COC','linewidth',1.5); 
    plot(rms(Referencia.Activo.y(:,10:12)),pisos,'--p','DisplayName','Activo','linewidth',1.5);
    xlabel('$\ddot{x}_{RMS}$ [m/s$^2$]','Interpreter','Latex')
    ylabel('Piso')
    grid on;
    xlim([0 1.1])
    yticks(pisos)
    xlm=xlim;
    title(terremotos(index),'Position',[xlm(2)/2 xlm(1)+3.05]);
    index = index+1;
end
lh = legend('Orientation','horizontal');
lh.Layout.Tile = 'South';
% exportgraphics(gcf,['Figs/CasoI/RespuestaRMS_casoI.jpg'],"Resolution",1000)

%% Tablas maximos 
for i = [Centro, Kobe, Maule]
    Referencia = i;
    pisos = [1,2,3];
    Filas = {'x1[mm]','x2[mm]','x3[mm]','d1[mm]','d2[mm]','d3[mm]',...
        'ddx1[mm/s^2]','ddx2[mm/s^2]','ddx3[mm/s^2]','F [N]'};

    piso = [pisos';pisos';pisos';1];

    maximos_sc = [(max(abs(Referencia.sincontrol.y(:,[1:6,10:12]))).*1000)';max(abs(Referencia.sincontrol.F))];
    maximos_poff = [(max(abs(Referencia.passiveoff.y(:,[1:6,10:12]))).*1000)';max(abs(Referencia.passiveoff.F))];
    maximos_pon = [(max(abs(Referencia.passiveon.y(:,[1:6,10:12]))).*1000)';max(abs(Referencia.passiveon.F))];
    maximos_LQG = [(max(abs(Referencia.LQGCOC.y(:,[1:6,10:12]))).*1000)';max(abs(Referencia.LQGCOC.F))];
    maximos_act = [(max(abs(Referencia.Activo.y(:,[1:6,10:12]))).*1000)';max(abs(Referencia.Activo.F))];
    energia_poff = 0;
    energia_pon = sum(Referencia.passiveon.ic)./2^12*12*4;
    energia_LQG = sum(Referencia.LQGCOC.ic)./2^12*12*4;
    
    Columnas = {'Piso','Sin control','Passive-Off','Passive-On','LQG-COC','Activo'};
    T = table(piso,maximos_sc,maximos_poff,maximos_pon,maximos_LQG,maximos_act,...
        'VariableNames',Columnas,'RowNames',Filas);

    % Porcentaje
    Filas = {'x1[%]','x2[%]','x3[%]','d1[%]','d2[%]','d3[%]',...
        'ddx1[%]','ddx2[%]','ddx3[%]','F [N]'};

    piso = [pisos';pisos';pisos';1];

    por_poff = (maximos_sc-maximos_poff)./maximos_sc*100;
    por_pon = (maximos_sc-maximos_pon)./maximos_sc*100;
    por_LQG = (maximos_sc-maximos_LQG)./maximos_sc*100;
    por_act = (maximos_sc-maximos_act)./maximos_sc*100;

    Columnas = {'Piso','Passive-Off','Passive-On','LQG-COC','Activo'};
    T = table(piso,por_poff,por_pon,por_LQG,por_act,...
    'VariableNames',Columnas,'RowNames',Filas)
end

%% Tablas RMS 
for i = [Centro, Kobe, Maule]
    Referencia = i;
    pisos = [1,2,3];
    Filas = {'x1[mm]','x2[mm]','x3[mm]','d1[mm]','d2[mm]','d3[mm]',...
        'ddx1[mm/s^2]','ddx2[mm/s^2]','ddx3[mm/s^2]','F [N]'};

    piso = [pisos';pisos';pisos';1];

    maximos_sc = [(rms(Referencia.sincontrol.y(:,[1:6,10:12])).*1000)';(rms(Referencia.sincontrol.F))];
    maximos_poff = [(rms(Referencia.passiveoff.y(:,[1:6,10:12])).*1000)';(rms(Referencia.passiveoff.F))];
    maximos_pon = [(rms(Referencia.passiveon.y(:,[1:6,10:12])).*1000)';(rms(Referencia.passiveon.F))];
    maximos_LQG = [(rms(Referencia.LQGCOC.y(:,[1:6,10:12])).*1000)';(rms(Referencia.LQGCOC.F))];
    maximos_act = [(rms(Referencia.Activo.y(:,[1:6,10:12])).*1000)';(rms(Referencia.Activo.F))];

    Columnas = {'Piso','Sin control','Passive-Off','Passive-On','LQG-COC','Activo'};
    T = table(piso,maximos_sc,maximos_poff,maximos_pon,maximos_LQG,maximos_act,...
        'VariableNames',Columnas,'RowNames',Filas);

    % Porcentaje
    Filas = {'x1[%]','x2[%]','x3[%]','d1[%]','d2[%]','d3[%]',...
        'ddx1[%]','ddx2[%]','ddx3[%]','F [N]'};

    piso = [pisos';pisos';pisos';1];

    por_poff = (maximos_sc-maximos_poff)./maximos_sc*100;
    por_pon = (maximos_sc-maximos_pon)./maximos_sc*100;
    por_LQG = (maximos_sc-maximos_LQG)./maximos_sc*100;
    por_act = (maximos_sc-maximos_act)./maximos_sc*100;

    Columnas = {'Piso','Passive-Off','Passive-On','LQG-COC','Activo'};
    T = table(piso,por_poff,por_pon,por_LQG,por_act,...
    'VariableNames',Columnas,'RowNames',Filas)
end

%% Grafico fuerza MR vs fuerza activa 
Referencia = load(['CasoI_Resultados/Centro/Centro50'],'Referencia').Referencia;

gcf = figure('Position', [10 10 1000 400]);
tiledlayout(2,1);
nexttile
plot(Referencia.t,Referencia.Activo.y(:,1).*1000,'k','DisplayName','Control activo ideal');
hold on;
plot(Referencia.t,Referencia.LQGCOC.y(:,1).*1000,'--r','DisplayName','Control semi-activo');
grid on;
xlabel('Tiempo [s]');
ylabel('Desp. [mm]');
xlim([10 30]);
lh = legend('Orientation','horizontal','Location','southeast');

nexttile
plot(Referencia.t,Referencia.Activo.F./1000,'k','DisplayName','Control activo ideal');
hold on;
plot(Referencia.t,Referencia.LQGCOC.F./1000,'--r','DisplayName','Control semi-activo');
grid on;
xlabel('Tiempo [s]');
ylabel('Fuerza [kN]');
xlim([10 30]);
% exportgraphics(gcf,['Figs/CasoI/Respuesta_casoI_FvsFA.jpg'],"Resolution",1000)

%% Gráfico Fuerza
gcf = figure('Position', [10 10 1000 500]);
plot(Referencia.t,Referencia.LQGCOC.F,'k','DisplayName','Clipped-Optimal Control');
hold on;
plot(Referencia.t,Referencia.passiveon.F,'r','DisplayName','Passive-On');
plot(Referencia.t,Referencia.passiveoff.F,'b','DisplayName','Passive-Off');
grid on;
xlabel('Tiempo [s]');
ylabel('Fuerza [N]');
xlim([0 Referencia.t(end)]);
legend();

%% Gráfico sistema sin control vs semi-activo
Referencia = load(['CasoI_Resultados/Centro/Centro50'],'Referencia').Referencia;

x3sc = Referencia.sincontrol.y(:,3);
x3cc = Referencia.LQGCOC.y(:,3);
x3max = (max(abs(x3sc))-max(abs(x3cc)))./max(abs(x3sc))*100;
x3rms = (rms(x3sc)-rms(x3cc))./rms(x3sc)*100;
d3sc = Referencia.sincontrol.y(:,6);
d3cc = Referencia.LQGCOC.y(:,6);
d3max = (max(abs(d3sc))-max(abs(d3cc)))./max(abs(d3sc))*100;
d3rms = (rms(d3sc)-rms(d3cc))./rms(d3sc)*100;
ddx3sc = Referencia.sincontrol.y(:,12);
ddx3cc = Referencia.LQGCOC.y(:,12);
ddx3max = (max(abs(ddx3sc))-max(abs(ddx3cc)))./max(abs(ddx3sc))*100;
ddx3rms = (rms(ddx3sc)-rms(ddx3cc))./rms(ddx3sc)*100;
t = Referencia.t;

gcf = figure('Position', [10 10 1000 400]);
tiledlayout(2,1);
hf = nexttile;
plot(t,x3sc.*1000,'k','DisplayName','Sin control');
hold on;
plot(t,x3cc.*1000,'--r','DisplayName','Control semi-activo');
grid on;
xlabel('Tiempo [s]');
ylabel('$x_3$ [mm]','Interpreter','Latex');
xlim([0 Referencia.t(end)]);
text(hf.XLim(end)*0.7,hf.YLim(end)*0.85,['Reducción del maximo de ',num2str(x3max,2),' %'],'FontSize',8)
text(hf.XLim(end)*0.7,hf.YLim(end)*0.65,['Reducción del RMS de ',num2str(x3rms,2),' %'],'FontSize',8)
lh = legend('Orientation','horizontal','Location','southeast');

hf = nexttile;
plot(t,ddx3sc.*1000,'k','DisplayName','Sin control');
hold on;
plot(t,ddx3cc.*1000,'--r','DisplayName','Control semi-activo');
grid on;
xlabel('Tiempo [s]');
ylabel('$\ddot{x}_3$ [mm/s$^2$]','Interpreter','Latex');
xlim([0 Referencia.t(end)]);
text(hf.XLim(end)*0.7,hf.YLim(end)*0.85,['Reducción del maximo de ',num2str(ddx3max,2),' %'],'FontSize',8)
text(hf.XLim(end)*0.7,hf.YLim(end)*0.65,['Reducción del RMS de ',num2str(ddx3rms,2),' %'],'FontSize',8)
% exportgraphics(gcf,['Figs/CasoI/Respuesta_casoI_RscvsRcc.jpg'],"Resolution",1000)


%% Cargar datos Centro diferentes magnitudes
% Escala 10%
Referencia = load(['CasoI_Resultados/Centro/Centro10'],'Referencia').Referencia;

maximos_sc = [max(abs(Referencia.sincontrol.y(:,[1:3,5:6,10:12])).*1000)';max(abs(Referencia.sincontrol.F))];
maximos_pon = [max(abs(Referencia.passiveon.y(:,[1:3,5:6,10:12])).*1000)';max(abs(Referencia.passiveon.F));sum(Referencia.passiveon.ic)*12*4/2^12];
maximos_LQG = [max(abs(Referencia.LQGCOC.y(:,[1:3,5:6,10:12])).*1000)';max(abs(Referencia.LQGCOC.F));sum(Referencia.LQGCOC.ic)*12*4/2^12];

Centro10_pon = [((maximos_sc(1:8)-maximos_pon(1:8))./maximos_sc(1:8)*100);maximos_pon(9:10)];
Centro10_LQG = [((maximos_sc(1:8)-maximos_LQG(1:8))./maximos_sc(1:8)*100);maximos_LQG(9:10)];

% Escala 30%
Referencia = load(['CasoI_Resultados/Centro/Centro30'],'Referencia').Referencia;

maximos_sc = [max(abs(Referencia.sincontrol.y(:,[1:3,5:6,10:12])).*1000)';max(abs(Referencia.sincontrol.F))];
maximos_pon = [max(abs(Referencia.passiveon.y(:,[1:3,5:6,10:12])).*1000)';max(abs(Referencia.passiveon.F));sum(Referencia.passiveon.ic)*12*4/2^12];
maximos_LQG = [max(abs(Referencia.LQGCOC.y(:,[1:3,5:6,10:12])).*1000)';max(abs(Referencia.LQGCOC.F));sum(Referencia.LQGCOC.ic)*12*4/2^12];

Centro30_pon = [((maximos_sc(1:8)-maximos_pon(1:8))./maximos_sc(1:8)*100);maximos_pon(9:10)];
Centro30_LQG = [((maximos_sc(1:8)-maximos_LQG(1:8))./maximos_sc(1:8)*100);maximos_LQG(9:10)];

% Escala 60%
Referencia = load(['CasoI_Resultados/Centro/Centro60'],'Referencia').Referencia;

maximos_sc = [max(abs(Referencia.sincontrol.y(:,[1:3,5:6,10:12])).*1000)';max(abs(Referencia.sincontrol.F))];
maximos_pon = [max(abs(Referencia.passiveon.y(:,[1:3,5:6,10:12])).*1000)';max(abs(Referencia.passiveon.F));sum(Referencia.passiveon.ic)*12*4/2^12];
maximos_LQG = [max(abs(Referencia.LQGCOC.y(:,[1:3,5:6,10:12])).*1000)';max(abs(Referencia.LQGCOC.F));sum(Referencia.LQGCOC.ic)*12*4/2^12];

Centro60_pon = [((maximos_sc(1:8)-maximos_pon(1:8))./maximos_sc(1:8)*100);maximos_pon(9:10)];
Centro60_LQG = [((maximos_sc(1:8)-maximos_LQG(1:8))./maximos_sc(1:8)*100);maximos_LQG(9:10)];

%% Cargar datos Kobe diferentes magnitudes
% Escala 10%
Referencia = load(['CasoI_Resultados/Kobe/Kobe10'],'Referencia').Referencia;

maximos_sc = [max(abs(Referencia.sincontrol.y(:,[1:3,5:6,10:12])).*1000)';max(abs(Referencia.sincontrol.F))];
maximos_pon = [max(abs(Referencia.passiveon.y(:,[1:3,5:6,10:12])).*1000)';max(abs(Referencia.passiveon.F));sum(Referencia.passiveon.ic)*12*4/2^12];
maximos_LQG = [max(abs(Referencia.LQGCOC.y(:,[1:3,5:6,10:12])).*1000)';max(abs(Referencia.LQGCOC.F));sum(Referencia.LQGCOC.ic)*12*4/2^12];

Kobe10_pon = [((maximos_sc(1:8)-maximos_pon(1:8))./maximos_sc(1:8)*100);maximos_pon(9:10)];
Kobe10_LQG = [((maximos_sc(1:8)-maximos_LQG(1:8))./maximos_sc(1:8)*100);maximos_LQG(9:10)];

% Escala 30%
Referencia = load(['CasoI_Resultados/Kobe/Kobe30'],'Referencia').Referencia;

maximos_sc = [max(abs(Referencia.sincontrol.y(:,[1:3,5:6,10:12])).*1000)';max(abs(Referencia.sincontrol.F))];
maximos_pon = [max(abs(Referencia.passiveon.y(:,[1:3,5:6,10:12])).*1000)';max(abs(Referencia.passiveon.F));sum(Referencia.passiveon.ic)*12*4/2^12];
maximos_LQG = [max(abs(Referencia.LQGCOC.y(:,[1:3,5:6,10:12])).*1000)';max(abs(Referencia.LQGCOC.F));sum(Referencia.LQGCOC.ic)*12*4/2^12];

Kobe30_pon = [((maximos_sc(1:8)-maximos_pon(1:8))./maximos_sc(1:8)*100);maximos_pon(9:10)];
Kobe30_LQG = [((maximos_sc(1:8)-maximos_LQG(1:8))./maximos_sc(1:8)*100);maximos_LQG(9:10)];

% Escala 60%
Referencia = load(['CasoI_Resultados/Kobe/Kobe60'],'Referencia').Referencia;

maximos_sc = [max(abs(Referencia.sincontrol.y(:,[1:3,5:6,10:12])).*1000)';max(abs(Referencia.sincontrol.F))];
maximos_pon = [max(abs(Referencia.passiveon.y(:,[1:3,5:6,10:12])).*1000)';max(abs(Referencia.passiveon.F));sum(Referencia.passiveon.ic)*12*4/2^12];
maximos_LQG = [max(abs(Referencia.LQGCOC.y(:,[1:3,5:6,10:12])).*1000)';max(abs(Referencia.LQGCOC.F));sum(Referencia.LQGCOC.ic)*12*4/2^12];

Kobe60_pon = [((maximos_sc(1:8)-maximos_pon(1:8))./maximos_sc(1:8)*100);maximos_pon(9:10)];
Kobe60_LQG = [((maximos_sc(1:8)-maximos_LQG(1:8))./maximos_sc(1:8)*100);maximos_LQG(9:10)];

%% Cargar datos Maule diferentes magnitudes
% Escala 10%
Referencia = load(['CasoI_Resultados/Maule/Maule10'],'Referencia').Referencia;

maximos_sc = [max(abs(Referencia.sincontrol.y(:,[1:3,5:6,10:12])).*1000)';max(abs(Referencia.sincontrol.F))];
maximos_pon = [max(abs(Referencia.passiveon.y(:,[1:3,5:6,10:12])).*1000)';max(abs(Referencia.passiveon.F));sum(Referencia.passiveon.ic)*12*4/2^12];
maximos_LQG = [max(abs(Referencia.LQGCOC.y(:,[1:3,5:6,10:12])).*1000)';max(abs(Referencia.LQGCOC.F));sum(Referencia.LQGCOC.ic)*12*4/2^12];

Maule10_pon = [((maximos_sc(1:8)-maximos_pon(1:8))./maximos_sc(1:8)*100);maximos_pon(9:10)];
Maule10_LQG = [((maximos_sc(1:8)-maximos_LQG(1:8))./maximos_sc(1:8)*100);maximos_LQG(9:10)];

% Escala 30%
Referencia = load(['CasoI_Resultados/Maule/Maule30'],'Referencia').Referencia;

maximos_sc = [max(abs(Referencia.sincontrol.y(:,[1:3,5:6,10:12])).*1000)';max(abs(Referencia.sincontrol.F))];
maximos_pon = [max(abs(Referencia.passiveon.y(:,[1:3,5:6,10:12])).*1000)';max(abs(Referencia.passiveon.F));sum(Referencia.passiveon.ic)*12*4/2^12];
maximos_LQG = [max(abs(Referencia.LQGCOC.y(:,[1:3,5:6,10:12])).*1000)';max(abs(Referencia.LQGCOC.F));sum(Referencia.LQGCOC.ic)*12*4/2^12];

Maule30_pon = [((maximos_sc(1:8)-maximos_pon(1:8))./maximos_sc(1:8)*100);maximos_pon(9:10)];
Maule30_LQG = [((maximos_sc(1:8)-maximos_LQG(1:8))./maximos_sc(1:8)*100);maximos_LQG(9:10)];

% Escala 60%
Referencia = load(['CasoI_Resultados/Maule/Maule60'],'Referencia').Referencia;

maximos_sc = [max(abs(Referencia.sincontrol.y(:,[1:3,5:6,10:12])).*1000)';max(abs(Referencia.sincontrol.F))];
maximos_pon = [max(abs(Referencia.passiveon.y(:,[1:3,5:6,10:12])).*1000)';max(abs(Referencia.passiveon.F));sum(Referencia.passiveon.ic)*12*4/2^12];
maximos_LQG = [max(abs(Referencia.LQGCOC.y(:,[1:3,5:6,10:12])).*1000)';max(abs(Referencia.LQGCOC.F));sum(Referencia.LQGCOC.ic)*12*4/2^12];

Maule60_pon = [((maximos_sc(1:8)-maximos_pon(1:8))./maximos_sc(1:8)*100);maximos_pon(9:10)];
Maule60_LQG = [((maximos_sc(1:8)-maximos_LQG(1:8))./maximos_sc(1:8)*100);maximos_LQG(9:10)];

%% Tabla diferentes magnitudes

Filas = {'El Centro 10%, passive-on','El Centro 10%, COC','El Centro 30%, passive-on','El Centro 30%, COC','El Centro 60%, passive-on','El Centro 60%, COC',...
    'Kobe 10%, passive-on','Kobe 10%, COC','Kobe 30%, passive-on','Kobe 30%, COC','Kobe 60%, passive-on','Kobe 60%, COC',...
    'Maule 10%, passive-on','Maule 10%, COC','Maule 30%, passive-on','Maule 30%, COC','Maule 60%, passive-on','Maule 60%, COC'};
Columnas = {'x_1','x_2','x_3','x_2-x_1','x_3-x_2','ddx_1','ddx_2','ddx_3','F[N]','E[VAs]'};

aux = [Centro10_pon,Centro10_LQG,Centro30_pon,Centro30_LQG,Centro60_pon,Centro60_LQG,...
    Kobe10_pon,Kobe10_LQG,Kobe30_pon,Kobe30_LQG,Kobe60_pon,Kobe60_LQG,...
    Maule10_pon,Maule10_LQG,Maule30_pon,Maule30_LQG,Maule60_pon,Maule60_LQG];
aux = aux';
T = array2table(aux,'VariableNames',Columnas,'RowNames',Filas)

%% Comparacion 3 casos
CasoI_centro = load(['CasoI_Resultados/Centro/Centro50'],'Referencia').Referencia;
CasoI_kobe = load(['CasoI_Resultados/Kobe/Kobe50'],'Referencia').Referencia;
CasoI_maule = load(['CasoI_Resultados/Maule/Maule50'],'Referencia').Referencia;

CasoII_centro = load(['CasoII_Resultados/Centro/Centro50'],'Referencia').Referencia;
CasoII_kobe = load(['CasoII_Resultados/Kobe/Kobe50'],'Referencia').Referencia;
CasoII_maule = load(['CasoII_Resultados/Maule/Maule50'],'Referencia').Referencia;

CasoIII_centro = load(['CasoIII_Resultados/Centro/Centro50'],'Referencia').Referencia;
CasoIII_kobe = load(['CasoIII_Resultados/Kobe/Kobe50'],'Referencia').Referencia;
CasoIII_maule = load(['CasoIII_Resultados/Maule/Maule50'],'Referencia').Referencia;

%% x3 max
x3_casoI = [max(abs(CasoI_centro.sincontrol.y(:,3))).*1000, max(abs(CasoI_centro.passiveoff.y(:,3))).*1000,max(abs(CasoI_centro.passiveon.y(:,3))).*1000,max(abs(CasoI_centro.LQGCOC.y(:,3))).*1000,max(abs(CasoI_centro.Activo.y(:,3))).*1000];
x3_casoII = [max(abs(CasoII_centro.sincontrol.y(:,3))).*1000, max(abs(CasoII_centro.passiveoff.y(:,3))).*1000,max(abs(CasoII_centro.passiveon.y(:,3))).*1000,max(abs(CasoII_centro.LQGCOC.y(:,3))).*1000,max(abs(CasoII_centro.Activo.y(:,3))).*1000];
x3_casoIII = [max(abs(CasoIII_centro.sincontrol.y(:,3))).*1000, max(abs(CasoIII_centro.passiveoff.y(:,3))).*1000,max(abs(CasoIII_centro.passiveon.y(:,3))).*1000,max(abs(CasoIII_centro.LQGCOC.y(:,3))).*1000,max(abs(CasoIII_centro.Activo.y(:,3))).*1000];

gcf = figure('Position', [10 10 1000 500]);
tiledlayout(2,3);
nexttile
bar([x3_casoI;x3_casoII;x3_casoIII]')
ylabel('$x_{3_{max}}$ [mm]','Interpreter','Latex','FontSize',14)
xticklabels({'Sin control','Passive-off','Passive-on','Semi-activo','Activo ideal'})
legend('Caso I','Caso II','Caso III');
grid on;
ylim([0 20])
title('El Centro 50%')

x3_casoI = [max(abs(CasoI_kobe.sincontrol.y(:,3))).*1000, max(abs(CasoI_kobe.passiveoff.y(:,3))).*1000,max(abs(CasoI_kobe.passiveon.y(:,3))).*1000,max(abs(CasoI_kobe.LQGCOC.y(:,3))).*1000,max(abs(CasoI_kobe.Activo.y(:,3))).*1000];
x3_casoII = [max(abs(CasoII_kobe.sincontrol.y(:,3))).*1000, max(abs(CasoII_kobe.passiveoff.y(:,3))).*1000,max(abs(CasoII_kobe.passiveon.y(:,3))).*1000,max(abs(CasoII_kobe.LQGCOC.y(:,3))).*1000,max(abs(CasoII_kobe.Activo.y(:,3))).*1000];
x3_casoIII = [max(abs(CasoIII_kobe.sincontrol.y(:,3))).*1000, max(abs(CasoIII_kobe.passiveoff.y(:,3))).*1000,max(abs(CasoIII_kobe.passiveon.y(:,3))).*1000,max(abs(CasoIII_kobe.LQGCOC.y(:,3))).*1000,max(abs(CasoIII_kobe.Activo.y(:,3))).*1000];
nexttile
bar([x3_casoI;x3_casoII;x3_casoIII]')
ylabel('$x_{3_{max}}$ [mm]','Interpreter','Latex','FontSize',14)
xticklabels({'Sin control','Passive-off','Passive-on','Semi-activo','Activo ideal'})
% legend('Caso I','Caso II','Caso III');
grid on;
ylim([0 20])
title('Kobe 50%')

x3_casoI = [max(abs(CasoI_maule.sincontrol.y(:,3))).*1000, max(abs(CasoI_maule.passiveoff.y(:,3))).*1000,max(abs(CasoI_maule.passiveon.y(:,3))).*1000,max(abs(CasoI_maule.LQGCOC.y(:,3))).*1000,max(abs(CasoI_maule.Activo.y(:,3))).*1000];
x3_casoII = [max(abs(CasoII_maule.sincontrol.y(:,3))).*1000, max(abs(CasoII_maule.passiveoff.y(:,3))).*1000,max(abs(CasoII_maule.passiveon.y(:,3))).*1000,max(abs(CasoII_maule.LQGCOC.y(:,3))).*1000,max(abs(CasoII_maule.Activo.y(:,3))).*1000];
x3_casoIII = [max(abs(CasoIII_maule.sincontrol.y(:,3))).*1000, max(abs(CasoIII_maule.passiveoff.y(:,3))).*1000,max(abs(CasoIII_maule.passiveon.y(:,3))).*1000,max(abs(CasoIII_maule.LQGCOC.y(:,3))).*1000,max(abs(CasoIII_maule.Activo.y(:,3))).*1000];
nexttile
bar([x3_casoI;x3_casoII;x3_casoIII]')
ylabel('$x_{3_{max}}$ [mm]','Interpreter','Latex','FontSize',14)
xticklabels({'Sin control','Passive-off','Passive-on','Semi-activo','Activo ideal'})
% legend('Caso I','Caso II','Caso III');
grid on;
ylim([0 20])
title('Maule 50%')

%% ddx3 max
x3_casoI = [max(abs(CasoI_centro.sincontrol.y(:,12))).*1000, max(abs(CasoI_centro.passiveoff.y(:,12))).*1000,max(abs(CasoI_centro.passiveon.y(:,12))).*1000,max(abs(CasoI_centro.LQGCOC.y(:,12))).*1000,max(abs(CasoI_centro.Activo.y(:,12))).*1000];
x3_casoII = [max(abs(CasoII_centro.sincontrol.y(:,12))).*1000, max(abs(CasoII_centro.passiveoff.y(:,12))).*1000,max(abs(CasoII_centro.passiveon.y(:,12))).*1000,max(abs(CasoII_centro.LQGCOC.y(:,12))).*1000,max(abs(CasoII_centro.Activo.y(:,12))).*1000];
x3_casoIII = [max(abs(CasoIII_centro.sincontrol.y(:,12))).*1000, max(abs(CasoIII_centro.passiveoff.y(:,12))).*1000,max(abs(CasoIII_centro.passiveon.y(:,12))).*1000,max(abs(CasoIII_centro.LQGCOC.y(:,12))).*1000,max(abs(CasoIII_centro.Activo.y(:,12))).*1000];

nexttile
bar([x3_casoI;x3_casoII;x3_casoIII]')
ylabel('$\ddot{x}_{3_{max}}$ [mm/s$^2$]','Interpreter','Latex','FontSize',14)
xticklabels({'Sin control','Passive-off','Passive-on','Semi-activo','Activo ideal'})
% legend('Caso I','Caso II','Caso III');
grid on;
ylim([0 8500])
title('El Centro 50%')

x3_casoI = [max(abs(CasoI_kobe.sincontrol.y(:,12))).*1000, max(abs(CasoI_kobe.passiveoff.y(:,12))).*1000,max(abs(CasoI_kobe.passiveon.y(:,12))).*1000,max(abs(CasoI_kobe.LQGCOC.y(:,12))).*1000,max(abs(CasoI_kobe.Activo.y(:,12))).*1000];
x3_casoII = [max(abs(CasoII_kobe.sincontrol.y(:,12))).*1000, max(abs(CasoII_kobe.passiveoff.y(:,12))).*1000,max(abs(CasoII_kobe.passiveon.y(:,12))).*1000,max(abs(CasoII_kobe.LQGCOC.y(:,12))).*1000,max(abs(CasoII_kobe.Activo.y(:,12))).*1000];
x3_casoIII = [max(abs(CasoIII_kobe.sincontrol.y(:,12))).*1000, max(abs(CasoIII_kobe.passiveoff.y(:,12))).*1000,max(abs(CasoIII_kobe.passiveon.y(:,12))).*1000,max(abs(CasoIII_kobe.LQGCOC.y(:,12))).*1000,max(abs(CasoIII_kobe.Activo.y(:,12))).*1000];
nexttile
bar([x3_casoI;x3_casoII;x3_casoIII]')
ylabel('$\ddot{x}_{3_{max}}$ [mm/s$^2$]','Interpreter','Latex','FontSize',14)
xticklabels({'Sin control','Passive-off','Passive-on','Semi-activo','Activo ideal'})
% legend('Caso I','Caso II','Caso III');
grid on;
ylim([0 8500])
title('Kobe 50%')

x3_casoI = [max(abs(CasoI_maule.sincontrol.y(:,12))).*1000, max(abs(CasoI_maule.passiveoff.y(:,12))).*1000,max(abs(CasoI_maule.passiveon.y(:,12))).*1000,max(abs(CasoI_maule.LQGCOC.y(:,12))).*1000,max(abs(CasoI_maule.Activo.y(:,12))).*1000];
x3_casoII = [max(abs(CasoII_maule.sincontrol.y(:,12))).*1000, max(abs(CasoII_maule.passiveoff.y(:,12))).*1000,max(abs(CasoII_maule.passiveon.y(:,12))).*1000,max(abs(CasoII_maule.LQGCOC.y(:,12))).*1000,max(abs(CasoII_maule.Activo.y(:,12))).*1000];
x3_casoIII = [max(abs(CasoIII_maule.sincontrol.y(:,12))).*1000, max(abs(CasoIII_maule.passiveoff.y(:,12))).*1000,max(abs(CasoIII_maule.passiveon.y(:,12))).*1000,max(abs(CasoIII_maule.LQGCOC.y(:,12))).*1000,max(abs(CasoIII_maule.Activo.y(:,12))).*1000];
nexttile
bar([x3_casoI;x3_casoII;x3_casoIII]')
ylabel('$\ddot{x}_{3_{max}}$ [mm/s$^2$]','Interpreter','Latex','FontSize',14)
xticklabels({'Sin control','Passive-off','Passive-on','Semi-activo','Activo ideal'})
% legend('Caso I','Caso II','Caso III');
grid on;
ylim([0 8500])
title('Maule 50%')
% exportgraphics(gcf,['Figs/x3_3casos.jpg'],"Resolution",1000)
