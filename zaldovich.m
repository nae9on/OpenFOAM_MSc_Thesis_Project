function zaldovich()
clc
clear all
close all

A = 1.8*10^8;
Ea = 38370;

T = 1600:1:2500;
k = A*exp(-Ea./T); 

T_trunc = T(1:100:end);
k_trunc = k(1:100:end);

figure();
hold on
plot(T,k);
plot(T(1:100:end),k(1:100:end),'*r');
xlabel('Temperature');
ylabel('Forward Rate Constant kf1')

%strValues = strtrim(cellstr(num2str([T_trunc(:) k_trunc(:)],'(%f,%f)')));
%text(T(1:100:end),k(1:100:end),strValues,'VerticalAlignment','bottom');
end

