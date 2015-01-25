clc
close all

if(1)
    
    path = '../axiBFR/postProcessing/sets/10000/';

    axial03 = load(strcat(path,'axial03_Ux.xy'));
    axial15 = load(strcat(path,'axial15_Ux.xy'));
    axial20 = load(strcat(path,'axial20_Ux.xy'));
    axial50 = load(strcat(path,'axial50_Ux.xy'));


    figure(1);
    subplot(2,2,1)
    plot(axial03(:,1),axial03(:,2),'-.ob');
    xlabel('Radial coordinate in m');
    ylabel('axial Velocity in m/s')
    title('axial velocity at x = 0.03m');
    subplot(2,2,2)
    plot(axial15(:,1),axial15(:,2),'-.ob');
    xlabel('Radial coordinate in m');
    ylabel('axial Velocity in m/s')
    title('axial velocity at x = 0.15m');
    subplot(2,2,3)
    plot(axial20(:,1),axial20(:,2),'-.ob');
    xlabel('Radial coordinate in m');
    ylabel('axial Velocity in m/s')
    title('axial velocity at x = 0.20m');
    subplot(2,2,4)
    plot(axial50(:,1),axial50(:,2),'-.ob');
    xlabel('Radial coordinate in m');
    ylabel('axial Velocity in m/s')
    title('axial velocity at x = 0.50m');

    radial03 = load(strcat(path,'radial03_Uy.xy'));
    radial15 = load(strcat(path,'radial15_Uy.xy'));
    radial20 = load(strcat(path,'radial20_Uy.xy'));
    radial50 = load(strcat(path,'radial50_Uy.xy'));


    figure(2);
    subplot(2,2,1)
    plot(radial03(:,1),radial03(:,2),'-.ob');
    xlabel('Radial coordinate in m');
    ylabel('radial Velocity in m/s')
    title('radial velocity at x = 0.03m');
    subplot(2,2,2)
    plot(radial15(:,1),radial15(:,2),'-.ob');
    xlabel('Radial coordinate in m');
    ylabel('radial Velocity in m/s')
    title('radial velocity at x = 0.15m');
    subplot(2,2,3)
    plot(radial20(:,1),radial20(:,2),'-.ob');
    xlabel('Radial coordinate in m');
    ylabel('radial Velocity in m/s')
    title('radial velocity at x = 0.20m');
    subplot(2,2,4)
    plot(radial50(:,1),radial50(:,2),'-.ob');
    xlabel('Radial coordinate in m');
    ylabel('radial Velocity in m/s')
    title('radial velocity at x = 0.50m');

    tangential03 = load(strcat(path,'tangential03_Uz.xy'));
    tangential15 = load(strcat(path,'tangential15_Uz.xy'));
    tangential20 = load(strcat(path,'tangential20_Uz.xy'));
    tangential50 = load(strcat(path,'tangential50_Uz.xy'));


    figure(3);
    subplot(2,2,1)
    plot(tangential03(:,1),tangential03(:,2),'-.ob');
    xlabel('Radial coordinate in m');
    ylabel('tangential Velocity in m/s')
    title('tangential velocity at x = 0.03m');
    subplot(2,2,2)
    plot(tangential15(:,1),tangential15(:,2),'-.ob');
    xlabel('Radial coordinate in m');
    ylabel('tangential Velocity in m/s')
    title('tangential velocity at x = 0.15m');
    subplot(2,2,3)
    plot(tangential20(:,1),tangential20(:,2),'-.ob');
    xlabel('Radial coordinate in m');
    ylabel('tangential Velocity in m/s')
    title('tangential velocity at x = 0.20m');
    subplot(2,2,4)
    plot(tangential50(:,1),tangential50(:,2),'-.ob');
    xlabel('Radial coordinate in m');
    ylabel('tangential Velocity in m/s')
    title('tangential velocity at x = 0.50m');

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fluent Results%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if(1)
    
    path = '../fluentResults/';

    axial03 = load(strcat(path,'axial03'));
    axial15 = load(strcat(path,'axial15'));
    axial20 = load(strcat(path,'axial20'));
    axial50 = load(strcat(path,'axial50'));


    figure(1);
    subplot(2,2,1); hold on
    plot(axial03(:,1),axial03(:,2),'r');
    xlabel('Radial coordinate in m');
    ylabel('axial Velocity in m/s')
    title('axial velocity at x = 0.03m');
    subplot(2,2,2); hold on
    plot(axial15(:,1),axial15(:,2),'r');
    xlabel('Radial coordinate in m');
    ylabel('axial Velocity in m/s')
    title('axial velocity at x = 0.15m');
    subplot(2,2,3); hold on
    plot(axial20(:,1),axial20(:,2),'r');
    xlabel('Radial coordinate in m');
    ylabel('axial Velocity in m/s')
    title('axial velocity at x = 0.20m');
    subplot(2,2,4); hold on
    plot(axial50(:,1),axial50(:,2),'r');
    xlabel('Radial coordinate in m');
    ylabel('axial Velocity in m/s')
    title('axial velocity at x = 0.50m');
    
    
    radial03 = load(strcat(path,'radial03'));
    radial15 = load(strcat(path,'radial15'));
    radial20 = load(strcat(path,'radial20'));
    radial50 = load(strcat(path,'radial50'));


    figure(2);
    subplot(2,2,1); hold on
    plot(radial03(:,1),radial03(:,2),'r');
    xlabel('Radial coordinate in m');
    ylabel('radial Velocity in m/s')
    title('radial velocity at x = 0.03m');
    subplot(2,2,2); hold on
    plot(radial15(:,1),radial15(:,2),'r');
    xlabel('Radial coordinate in m');
    ylabel('radial Velocity in m/s')
    title('radial velocity at x = 0.15m');
    subplot(2,2,3); hold on
    plot(radial20(:,1),radial20(:,2),'r');
    xlabel('Radial coordinate in m');
    ylabel('radial Velocity in m/s')
    title('radial velocity at x = 0.20m');
    subplot(2,2,4); hold on
    plot(radial50(:,1),radial50(:,2),'r');
    xlabel('Radial coordinate in m');
    ylabel('radial Velocity in m/s')
    title('radial velocity at x = 0.50m');

    tangential03 = load(strcat(path,'tangential03'));
    tangential15 = load(strcat(path,'tangential15'));
    tangential20 = load(strcat(path,'tangential20'));
    tangential50 = load(strcat(path,'tangential50'));


    figure(3);
    subplot(2,2,1); hold on
    plot(tangential03(:,1),-tangential03(:,2),'r');
    xlabel('Radial coordinate in m');
    ylabel('tangential Velocity in m/s')
    title('tangential velocity at x = 0.03m');
    subplot(2,2,2); hold on
    plot(tangential15(:,1),-tangential15(:,2),'r');
    xlabel('Radial coordinate in m');
    ylabel('tangential Velocity in m/s')
    title('tangential velocity at x = 0.15m');
    subplot(2,2,3); hold on
    plot(tangential20(:,1),-tangential20(:,2),'r');
    xlabel('Radial coordinate in m');
    ylabel('tangential Velocity in m/s')
    title('tangential velocity at x = 0.20m');
    subplot(2,2,4); hold on
    plot(tangential50(:,1),-tangential50(:,2),'r');
    xlabel('Radial coordinate in m');
    ylabel('tangential Velocity in m/s')
    title('tangential velocity at x = 0.50m');

end

figure(1); hold on
subplot(2,2,1);
legend('OpenFOAM','Fluent');
subplot(2,2,2);
legend('OpenFOAM','Fluent');
subplot(2,2,3);
legend('OpenFOAM','Fluent');
subplot(2,2,4);
legend('OpenFOAM','Fluent');

figure(2); hold on
subplot(2,2,1);
legend('OpenFOAM','Fluent');
subplot(2,2,2);
legend('OpenFOAM','Fluent');
subplot(2,2,3);
legend('OpenFOAM','Fluent');
subplot(2,2,4);
legend('OpenFOAM','Fluent');

figure(3); hold on
subplot(2,2,1);
legend('OpenFOAM','Fluent');
subplot(2,2,2);
legend('OpenFOAM','Fluent');
subplot(2,2,3);
legend('OpenFOAM','Fluent');
subplot(2,2,4);
legend('OpenFOAM','Fluent');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% path = '../3dBFR/postProcessing/sets/10000/';
% 
% axial03 = load(strcat(path,'axial03_Ux.xy'));
% axial15 = load(strcat(path,'axial15_Ux.xy'));
% axial20 = load(strcat(path,'axial20_Ux.xy'));
% axial50 = load(strcat(path,'axial50_Ux.xy'));
% 
% 
% figure(1);
% subplot(2,2,1); hold on
% plot(axial03(:,1),axial03(:,2),'r');
% xlabel('Radial coordinate in m');
% ylabel('axial Velocity in m/s')
% title('axial velocity at x = 0.03m');
% subplot(2,2,2); hold on
% plot(axial15(:,1),axial15(:,2),'r');
% xlabel('Radial coordinate in m');
% ylabel('axial Velocity in m/s')
% title('axial velocity at x = 0.15m');
% subplot(2,2,3); hold on
% plot(axial20(:,1),axial20(:,2),'r');
% xlabel('Radial coordinate in m');
% ylabel('axial Velocity in m/s')
% title('axial velocity at x = 0.20m');
% subplot(2,2,4); hold on
% plot(axial50(:,1),axial50(:,2),'r');
% xlabel('Radial coordinate in m');
% ylabel('axial Velocity in m/s')
% title('axial velocity at x = 0.50m');
% 
% radial03 = load(strcat(path,'radial03_Uy.xy'));
% radial15 = load(strcat(path,'radial15_Uy.xy'));
% radial20 = load(strcat(path,'radial20_Uy.xy'));
% radial50 = load(strcat(path,'radial50_Uy.xy'));
% 
% 
% figure(2);
% subplot(2,2,1); hold on
% plot(radial03(:,1),radial03(:,2),'r');
% xlabel('Radial coordinate in m');
% ylabel('radial Velocity in m/s')
% title('radial velocity at x = 0.03m');
% subplot(2,2,2); hold on
% plot(radial15(:,1),radial15(:,2),'r');
% xlabel('Radial coordinate in m');
% ylabel('radial Velocity in m/s')
% title('radial velocity at x = 0.15m');
% subplot(2,2,3); hold on
% plot(radial20(:,1),radial20(:,2),'r');
% xlabel('Radial coordinate in m');
% ylabel('radial Velocity in m/s')
% title('radial velocity at x = 0.20m');
% subplot(2,2,4); hold on
% plot(radial50(:,1),radial50(:,2),'r');
% xlabel('Radial coordinate in m');
% ylabel('radial Velocity in m/s')
% title('radial velocity at x = 0.50m');
% 
% tangential03 = load(strcat(path,'tangential03_Uz.xy'));
% tangential15 = load(strcat(path,'tangential15_Uz.xy'));
% tangential20 = load(strcat(path,'tangential20_Uz.xy'));
% tangential50 = load(strcat(path,'tangential50_Uz.xy'));
% 
% 
% figure(3);
% subplot(2,2,1); hold on
% plot(tangential03(:,1),tangential03(:,2),'r');
% xlabel('Radial coordinate in m');
% ylabel('tangential Velocity in m/s')
% title('tangential velocity at x = 0.03m');
% subplot(2,2,2); hold on
% plot(tangential15(:,1),tangential15(:,2),'r');
% xlabel('Radial coordinate in m');
% ylabel('tangential Velocity in m/s')
% title('tangential velocity at x = 0.15m');
% subplot(2,2,3); hold on
% plot(tangential20(:,1),tangential20(:,2),'r');
% xlabel('Radial coordinate in m');
% ylabel('tangential Velocity in m/s')
% title('tangential velocity at x = 0.20m');
% subplot(2,2,4); hold on
% plot(tangential50(:,1),tangential50(:,2),'r');
% xlabel('Radial coordinate in m');
% ylabel('tangential Velocity in m/s')
% title('tangential velocity at x = 0.50m');
