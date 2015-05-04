    clc
    close all
    
    if(1)
        
        path = '../fluentResults/EDM/';

        CO2 = load(strcat(path,'axis_CO2'));
        H2O = load(strcat(path,'axis_H2O'));
        CH4 = load(strcat(path,'axis_CH4'));
        N2 = load(strcat(path,'axis_N2'));
        O2 = load(strcat(path,'axis_O2'));
        T15 = load(strcat(path,'at15_T'));
        T20 = load(strcat(path,'at20_T'));
        T50 = load(strcat(path,'at50_T'));
        
      
        figure(1);
        hold on

        plot(CO2(1:3:end,1),CO2(1:3:end,2),'-.or','MarkerSize',5)
        plot(H2O(1:3:end,1),H2O(1:3:end,2),'-.sb','MarkerSize',5)
        plot(CH4(1:3:end,1),CH4(1:3:end,2),'-.vk','MarkerSize',5)
        %plot(N2(:,1),N2(:,2),'-.ok','MarkerSize',5)
        plot(O2(1:3:end,1),O2(1:3:end,2),'-.+m','MarkerSize',5)

        xlabel('Axial coordinate (m)','FontSize', 30,'Color','k');
        ylabel('Mass Fraction','FontSize', 30,'Color','k')
        title('Mass Fraction of different species','FontSize', 30);
        h_legend = legend('CO2','H2O','CH4','O2');
        set(h_legend,'FontSize',30,'fontweight','bold');

        figure(2);
        hold on

        plot(T15(1:3:end,1),T15(1:3:end,2),'-.or','MarkerSize',5)
        plot(T50(1:3:end,1),T50(1:3:end,2),'bs','MarkerSize',5)

        xlabel('Radial coordinate (m)','FontSize', 15,'Color','k');
        ylabel('Temperature (K)','FontSize', 15,'Color','k')
        h_legend = legend('0.2m from inlet','0.5m from inlet');
        set(h_legend,'FontSize',15,'fontweight','bold');
        
        disp('Fluent, mass fractions at the outlet');
        disp(['CO2 = ',num2str(CO2(end))])
        disp(['H2O = ',num2str(H2O(end))])
        disp(['CH4 = ',num2str(CH4(end))])
        %disp(['N2 = ',num2str(N2(end))])
        disp(['O2 = ',num2str(O2(end))])
        
    end
    
    if(1)
        
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
        
        figure(1);
        hold on
        plot(CO2(1:3:end,1),CO2(1:3:end,2),'-.r','LineWidth',2,'MarkerSize',1)
        plot(H2O(1:3:end,1),H2O(1:3:end,2),'-.b','LineWidth',2,'MarkerSize',1)
        plot(CH4(1:3:end,1),CH4(1:3:end,2),'-.k','LineWidth',2,'MarkerSize',1)
        %plot(N2(:,1),N2(:,2),'-.ok','LineWidth',2,'MarkerSize',1)
        plot(O2(1:3:end,1),O2(1:3:end,2),'-.m','LineWidth',2,'MarkerSize',1)
        axis([0 0.5 0 1.1])
        xlabel('Axial coordinate (m)','FontSize', 15,'Color','k');
        ylabel('Mass Fraction','FontSize', 15,'Color','k')
        title('Mass Fraction of different species','FontSize', 15);
        h_legend = legend('CO2','H2O','CH4','O2');
        set(h_legend,'FontSize',15,'fontweight','bold');

        figure(2);
        hold on

        plot(T15(:,1),T15(:,2),'-.or','LineWidth',2,'MarkerSize',1)
        %plot(T50(:,1),T50(:,2),'-.ob','LineWidth',3,'MarkerSize',1)
        plot(T500(:,1),T500(:,2),'-.ob','LineWidth',2,'MarkerSize',1)

        xlabel('Radial coordinate (m)','FontSize', 15,'Color','k');
        ylabel('Temperature (K)','FontSize', 15,'Color','k')
        h_legend = legend('0.2m from inlet','0.5m from inlet');
        set(h_legend,'FontSize',15,'fontweight','bold');
    
        
        disp('OpenFOAM, mass fractions at the outlet');
        disp(['CO2 = ',num2str(CO2(end))])
        disp(['H2O = ',num2str(H2O(end))])
        disp(['CH4 = ',num2str(CH4(end))])
        %disp(['N2 = ',num2str(N2(end))])
        disp(['O2 = ',num2str(O2(end))])
        
    end
    
    

    
    
    if(0)
        
        path = '../fluentResults/EDC/';

        CO2 = load(strcat(path,'axis_CO2'));
        H2O = load(strcat(path,'axis_H2O'));
        CH4 = load(strcat(path,'axis_CH4'));
        N2 = load(strcat(path,'axis_N2'));
        O2 = load(strcat(path,'axis_O2'));
        T15 = load(strcat(path,'axisx15_T'));
        T50 = load(strcat(path,'axisx50_T'));
        T200 = load(strcat(path,'axisx200_T'));
        T500 = load(strcat(path,'axisx500_T'));
        
        figure(1);
        hold on

        plot(CO2(:,1),CO2(:,2),'-.sr','MarkerSize',3)
        plot(H2O(:,1),H2O(:,2),'-.sb','MarkerSize',3)
        plot(CH4(:,1),CH4(:,2),'-.sg','MarkerSize',3)
        plot(N2(:,1),N2(:,2),'-.sk','MarkerSize',3)
        plot(O2(:,1),O2(:,2),'-.sy','MarkerSize',3)
        axis([0 1 0 1])

        xlabel('Axial coordinate (m)','FontSize', 15,'Color','b');
        ylabel('Mass Fraction','FontSize', 15,'Color','b')
        title('Mass Fraction of different species','FontSize', 15);
        h_legend = legend('CO2','H2O','CH4','N2','O2');
        set(h_legend,'FontSize',14,'fontweight','bold');

        figure(2);
        hold on

        plot(T15(:,1),T15(:,2),'-.sr','MarkerSize',3)
        %plot(T50(:,1),T50(:,2),'-.sb','MarkerSize',5)
        plot(T500(:,1),T500(:,2),'-.sb','MarkerSize',3)

        xlabel('Radial coordinate (m)','FontSize', 15,'Color','k');
        ylabel('Temperature (K)','FontSize', 15,'Color','k')
        h_legend = legend('close to inlet','far from inlet');
        set(h_legend,'FontSize',14,'fontweight','bold');
    
    end
    
    
