%% Correr simulaciones
terremoto = 'Centro'; 
escala = 50;
caso = 'CasoI';
vRTHS(terremoto,escala,caso);

%% Cargar resultados
load(['Resultados/',caso,'/',terremoto,num2str(escala),'SinControl'],'SinControl')
load(['Resultados/',caso,'/',terremoto,num2str(escala),'Referencia'],'Referencia')
load(['Resultados/',caso,'/',terremoto,num2str(escala),'SinCompensar'],'SinCompensar')
load(['Resultados/',caso,'/',terremoto,num2str(escala),'ConCompensar'],'ConCompensar')

%% Obtener indicadores
indicador_sc_poff = indicadores(SinControl,Referencia.passiveoff,SinCompensar.passiveoff)';
indicador_cc_poff = indicadores(SinControl,Referencia.passiveoff,ConCompensar.passiveoff)';

indicador_sc_pon = indicadores(SinControl,Referencia.passiveon,SinCompensar.passiveon)';
indicador_cc_pon = indicadores(SinControl,Referencia.passiveon,ConCompensar.passiveon)';

indicador_sc_COC = indicadores(SinControl,Referencia.LQGCOC,SinCompensar.LQGCOC)';
indicador_cc_COC = indicadores(SinControl,Referencia.LQGCOC,ConCompensar.LQGCOC)';

aux = [indicador_sc_poff,indicador_cc_poff,indicador_sc_pon,indicador_cc_pon,indicador_sc_COC,indicador_cc_COC];


%% Graficos de los indicadores de cada caso
for caso = ["CasoI","CasoII","CasoIII"]
    indicador_sc_poff = {};
    indicador_cc_poff = {};
    indicador_sc_pon = {};
    indicador_cc_pon = {};
    indicador_sc_COC = {};
    indicador_cc_COC = {};
    escala = 50;

    i = 1;
    for terremoto = ["Centro","Kobe","Maule"]
        load(strcat('Resultados/',caso,'/',terremoto,num2str(escala),'SinControl'),'SinControl')
        load(strcat('Resultados/',caso,'/',terremoto,num2str(escala),'Referencia'),'Referencia')
        load(strcat('Resultados/',caso,'/',terremoto,num2str(escala),'SinCompensar'),'SinCompensar')
        load(strcat('Resultados/',caso,'/',terremoto,num2str(escala),'ConCompensar'),'ConCompensar')

        indicador_sc_poff{i} = indicadores(SinControl,Referencia.passiveoff,SinCompensar.passiveoff)';
        indicador_cc_poff{i} = indicadores(SinControl,Referencia.passiveoff,ConCompensar.passiveoff)';

        indicador_sc_pon{i} = indicadores(SinControl,Referencia.passiveon,SinCompensar.passiveon)';
        indicador_cc_pon{i} = indicadores(SinControl,Referencia.passiveon,ConCompensar.passiveon)';

        indicador_sc_COC{i} = indicadores(SinControl,Referencia.LQGCOC,SinCompensar.LQGCOC)';
        indicador_cc_COC{i} = indicadores(SinControl,Referencia.LQGCOC,ConCompensar.LQGCOC)';

        i = i+1;
    end
    max_indicador_sc_poff = max(cell2mat(indicador_sc_poff),[],2);
    max_indicador_cc_poff = max(cell2mat(indicador_cc_poff),[],2);
    max_indicador_sc_pon = max(cell2mat(indicador_sc_pon),[],2);
    max_indicador_cc_pon = max(cell2mat(indicador_cc_pon),[],2);
    max_indicador_sc_COC = max(cell2mat(indicador_sc_COC),[],2);
    max_indicador_cc_COC = max(cell2mat(indicador_cc_COC),[],2);

    gcf = figure('Position', [10 10 1000 300]);
    tiledlayout(1,3);
    nexttile;
    bar([max_indicador_sc_poff,max_indicador_cc_poff])
    set(gca,'xticklabel',["J_1","J_2","J_3","J_4","J_5","J_6","J_7","J_8","J_9","J_{10}"])
    set(gca,'YScale','log')
    ylim([10^-1 10^3.3])
    grid on;
    title('Passive-off')

    nexttile;
    bar([max_indicador_sc_pon,max_indicador_cc_pon])
    set(gca,'xticklabel',["J_1","J_2","J_3","J_4","J_5","J_6","J_7","J_8","J_9","J_{10}"])
    set(gca,'YScale','log')
    ylim([10^-1 10^3.3])
    grid on;
    title('Passive-on')

    nexttile;
    bar([max_indicador_sc_COC,max_indicador_cc_COC])
    set(gca,'xticklabel',["J_1","J_2","J_3","J_4","J_5","J_6","J_7","J_8","J_9","J_{10}"])
    
    set(gca,'YScale','log')
    ylim([10^-1 10^3.3])
    grid on;
    title('Semi-activo')
    
    legend('Sin compensación','Con compensación')
    lh = legend('Orientation','horizontal');
    lh.Layout.Tile = 'South';
%     exportgraphics(gcf,strcat('Figs/',caso,'_indicadores.jpg'),"Resolution",1000)
end