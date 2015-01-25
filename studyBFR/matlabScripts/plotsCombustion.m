clc
close all

if(1)
    
    path = '../combustionAxi/postProcessing/sets/10/';
    
    CO2 = load(strcat(path,'axis_CO2.xy'));
    H2O = load(strcat(path,'axis_H2O.xy'));
    CH4 = load(strcat(path,'axis_CH4.xy'));
    N2 = load(strcat(path,'axis_N2.xy'));
    O2 = load(strcat(path,'axis_O2.xy'));
    
    
    CO2(end,2)
    H2O(end,2)
    CH4(end,2)
    N2(end,2)
    O2(end,2)
    
    
   
    figure(1);
    subplot(2,2,1);
    plot(CO2(:,1),CO2(:,2),'-.ob')
    xlabel('Axial coordinate in m');
    ylabel('CO2 mass fraction')
    title('CO2');
    
    subplot(2,2,2);
    plot(H2O(:,1),H2O(:,2),'-.ob')
    xlabel('Axial coordinate in m');
    ylabel('H2O mass fraction')
    title('H2O');
    
    subplot(2,2,3);
    plot(CH4(:,1),CH4(:,2),'-.ob')
    xlabel('Axial coordinate in m');
    ylabel('CH4 mass fraction')
    title('CH4');
    
    subplot(2,2,4);
    plot(N2(:,1),N2(:,2),'-.ob')
    xlabel('Axial coordinate in m');
    ylabel('N2 mass fraction')
    title('N2');
    
    figure();
    plot(O2(:,1),O2(:,2),'-.ob')
    xlabel('Axial coordinate in m');
    ylabel('O2 mass fraction')
    title('O2');
    

    
end
