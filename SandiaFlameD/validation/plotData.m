clc
clear all
close all

% foamCalc components U
% foamCalc mag U

%parameters
preInlet = 14.58;
diameter = 7.2e-03;

%formating and paths
format shortG %shortEng compact
pathExp = './experimental/statistics/';
extension = '.Yave';

%reading experimental data
D01 = importdata(strcat(pathExp,'D01',extension));
D02 = importdata(strcat(pathExp,'D02',extension));
D03 = importdata(strcat(pathExp,'D03',extension));
D075 = importdata(strcat(pathExp,'D075',extension));
D15 = importdata(strcat(pathExp,'D15',extension));
D30 = importdata(strcat(pathExp,'D30',extension));
D45 = importdata(strcat(pathExp,'D45',extension));
D60 = importdata(strcat(pathExp,'D60',extension));
D75 = importdata(strcat(pathExp,'D75',extension));

index = strcat(num2str([1:max(size(D01.textdata(4,:)))]'),'__');
strcat(index,D01.textdata(4,:)')

xd = [1 2 3 7.5 15 30 45 60 75];
N = max(size(xd));

dataFiles = {D01,D02,D03,D075,D15,D30,D45,D60,D75};

% ' 1__r/d'
% ' 2__F'
% ' 3__Frms'
% ' 4__T(K)'
% ' 5__Trms'
% ' 6__YO2'
% ' 7__YO2rms'
% ' 8__YN2'
% ' 9__YN2rms'
% '10__YH2'
% '11__YH2rms'
% '12__YH2O'
% '13__YH2Orms'
% '14__YCH4'
% '15__YCH4rms'
% '16__YCO'
% '17__YCOrms'
% '18__YCO2'
% '19__YCO2rms'
% '20__YOH'
% '21__YOHrms'
% '22__YNO'
% '23__YNOrms'
% '24__YCOLIF'
% '25__YCOrms'
% '26__TNDR'

%extracting some quantity Q
QtyExp = 4;

%axial_CH4_CO2_H2O_N2_O2_T_Ux_Uy_Uz_magU
%axial_CH4_CO_CO2_H2_H2O_N2_NO_O2_OH_T_Ux_Uy_Uz_magU.xy
QtyOF = 7;

dataQ = [];
for i=1:N
    valMat = dataFiles{i}.data;
    dataQ = [dataQ; xd(i)*ones(max(size(valMat(:,1))),1) valMat(:,1) valMat(:,QtyExp)];
end

%extracting quantity along axis
dataQaxis = [];
for i=1:max(size(dataQ))
    if dataQ(i,2) == 0
        dataQaxis = [dataQaxis; dataQ(i,:)];
    end
end

pathOF_NoRadGLB = '../refined_myFlameD_GLB/postProcessing/sets/0.3/';
Qty_NoRadGLB = load(strcat(pathOF_NoRadGLB,'axial_CH4_CO2_H2O_N2_O2_T_Ux_Uy_Uz_magU.xy'));

pathOF_RadGLB = '../refined_myFlameD_GLB/postProcessing/sets/0.6/';
Qty_RadGLB = load(strcat(pathOF_RadGLB,'axial_CH4_CO2_H2O_N2_O2_T_Ux_Uy_Uz_magU.xy'));

pathOF_RadGRI3 = '../myFlameD_GRI3/postProcessing/sets/0.2/';
Qty_RadGRI3 = load(strcat(pathOF_RadGRI3,'axial_CH4_CO2_H2O_N2_O2_T_Ux_Uy_Uz_magU.xy'));

figure(1);
hold on
coarse = 10;

plot(dataQaxis(:,1), dataQaxis(:,3),'ok','MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',10);
plot(Qty_NoRadGLB(1:coarse:end,1)/diameter-preInlet,Qty_NoRadGLB(1:coarse:end,QtyOF),'-.sb','MarkerSize',5);
plot(Qty_RadGLB(1:coarse:end,1)/diameter-preInlet,Qty_RadGLB(1:coarse:end,QtyOF),'-.vr','MarkerSize',5);
plot(Qty_RadGRI3(1:coarse:end,1)/diameter-preInlet,Qty_RadGRI3(1:coarse:end,QtyOF),'-.ok','LineWidth',1,'MarkerSize',5);

xlabel('x/d','FontSize', 20,'Color','k');
ylabel(D01.textdata(4,QtyExp),'FontSize', 20,'Color','k');
h_legend = legend('Exp','ke-GLB','ke-GLB-P1','ke-GRI3-P1');
%h_legend = legend('Exp','EDC-GRI3');
set(h_legend,'FontSize',12,'fontweight','bold');
title('Sandia Flame D','FontSize', 15,'Color','k');


%%
%Radial Plots
Qty_NoRadGLB = load(strcat(pathOF_NoRadGLB,'radial_CH4_CO2_H2O_N2_O2_T_Ux_Uy_Uz_magU.xy'));
Qty_RadGLB = load(strcat(pathOF_RadGLB,'radial_CH4_CO2_H2O_N2_O2_T_Ux_Uy_Uz_magU.xy'));
Qty_RadGRI3 = load(strcat(pathOF_RadGRI3,'radial_CH4_CO2_H2O_N2_O2_T_Ux_Uy_Uz_magU.xy'));

coarse = 10;

figure(2);
hold on
plot(D45.data(:,1), D45.data(:,QtyExp),'ok','MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',10);
plot(Qty_NoRadGLB(1:coarse:end,1)/diameter,Qty_NoRadGLB(1:coarse:end,QtyOF),'-.sb','MarkerSize',5);
plot(Qty_RadGLB(1:coarse:end,1)/diameter,Qty_RadGLB(1:coarse:end,QtyOF),'-.vr','MarkerSize',5);
plot(Qty_RadGRI3(1:coarse:end,1)/diameter,Qty_RadGRI3(1:coarse:end,QtyOF),'-.ok','LineWidth',1,'MarkerSize',5);

xlabel('r/d','FontSize', 20,'Color','k');
ylabel(D45.textdata(4,QtyExp),'FontSize', 20,'Color','k');
h_legend = legend('Exp','ke-GLB','ke-GLB-P1','ke-GRI3-P1');
%h_legend = legend('Exp','EDC-GRI3');
set(h_legend,'FontSize',12,'fontweight','bold');
title('Sandia Flame D','FontSize', 15,'Color','k');
