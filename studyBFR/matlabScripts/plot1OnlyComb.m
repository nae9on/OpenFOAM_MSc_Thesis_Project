    clc
    close all
    
    path = '../combustionAxi/postProcessing/sets/10/';

    CO2 = load(strcat(path,'axis_CO2.xy'));
    H2O = load(strcat(path,'axis_H2O.xy'));
    CH4 = load(strcat(path,'axis_CH4.xy'));
    N2 = load(strcat(path,'axis_N2.xy'));
    O2 = load(strcat(path,'axis_O2.xy'));
    T15 = load(strcat(path,'axisx15_T.xy'));
    T50 = load(strcat(path,'axisx50_T.xy'));
    T200 = load(strcat(path,'axisx200_T.xy'));
    T500 = load(strcat(path,'axisx500_T.xy'));
    
    
    figure();
    hold on
    
    plot(CO2(:,1),CO2(:,2),'-.or','LineWidth',3,'MarkerSize',1)
    plot(H2O(:,1),H2O(:,2),'-.ob','LineWidth',3,'MarkerSize',1)
    plot(CH4(:,1),CH4(:,2),'-.og','LineWidth',3,'MarkerSize',1)
    plot(N2(:,1),N2(:,2),'-.ok','LineWidth',3,'MarkerSize',1)
    plot(O2(:,1),O2(:,2),'-.oy','LineWidth',3,'MarkerSize',1)

    xlabel('Axial coordinate (m)','FontSize', 15,'Color','b');
    ylabel('Mass Fraction','FontSize', 15,'Color','b')
    title('Mass Fraction of different species','FontSize', 15);
    h_legend = legend('CO2','H2O','CH4','N2','O2');
    set(h_legend,'FontSize',14,'fontweight','bold');
    
    figure();
    hold on
    
    plot(T15(:,1),T15(:,2),'-.or','LineWidth',3,'MarkerSize',1)
    %plot(T50(:,1),T50(:,2),'-.ob','LineWidth',3,'MarkerSize',1)
    plot(T500(:,1),T500(:,2),'-.ob','LineWidth',3,'MarkerSize',1)

    xlabel('Radial coordinate (m)','FontSize', 15,'Color','k');
    ylabel('Temperature (K)','FontSize', 15,'Color','k')
    h_legend = legend('close to inlet','far from inlet');
    set(h_legend,'FontSize',14,'fontweight','bold');
    
