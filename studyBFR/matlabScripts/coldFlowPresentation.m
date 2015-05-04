    clc
    close all
    
    path = '../axiBFR/postProcessing/sets/10000/';
    axial15 = load(strcat(path,'axial15_Ux.xy'));
    
    figure();
    hold on
    plot(axial15(:,1),axial15(:,2),'-.ob','LineWidth',1,'MarkerSize',10);
    
    path = '../fluentResults/coldFlow/';
    axial15 = load(strcat(path,'axial15'));
    plot(axial15(:,1),axial15(:,2),'r','LineWidth',2,'MarkerSize',10);
    
    xlabel('Radial coordinate (m)','FontSize', 15);
    ylabel('axial Velocity (m/s)','FontSize', 15)
    title('Axial Velocity 15 cm Away from the Inlet','FontSize', 15);
    h_legend = legend('OpenFOAM','Fluent');
    set(h_legend,'FontSize',14,'fontweight','bold');