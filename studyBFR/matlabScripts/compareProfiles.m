clc
clear all
close all

path0 = '../axiBFRold/postProcessing/sets/6000/';
path1 = '../2dBFR/postProcessing/sets/6000/';
path2 = '../axiBFR/postProcessing/sets/6000/';
path3 = '../3dBFR/postProcessing/sets/6000/';

Axial15_0 = load(strcat(path0,'axial15_Ux.xy'));
Axial15_1 = load(strcat(path1,'axial15_Ux.xy'));
Axial15_2 = load(strcat(path2,'axial15_Ux.xy'));
Axial15_3 = load(strcat(path3,'axial15_Ux.xy'));


figure();
hold on
plot(Axial15_0(:,1),Axial15_0(:,2),'-.ok');
plot(Axial15_1(:,1),Axial15_1(:,2),'-.ob');
plot(Axial15_2(:,1),Axial15_2(:,2),'-.og');
plot(Axial15_3(:,1),Axial15_3(:,2),'-.or');

legend('axi old','2d','axi new with effecient mesh(half the previous control volumes)','3d');
xlabel('Radial coordinate in m');
ylabel('Axial Velocity in m/s')