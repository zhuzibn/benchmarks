clear all;clc;tic

FL_length=25e-9;FL_width=25*pi*1e-9;FL_thickness=1.4e-9; %free layer geometry [m]
%% magnetic parameters
% alpha, Ku, Ms obtained from spice file of 'manual_Purdue Nanoelectronics ...
%Research Laboratory Magnetic Tunnel Junction Model'
alpha=0.028;
P1=0.8;P2=0.3;%[P1:pin layer, P2: free layer]
kB=1.38e-23;%[J.K-1]
Ea=56*kB*300;%[J]
%Ku=8.16e5; %erg/cm3
Ms=700; %emu/cm3
%% calculation
FL_volume=FL_length*FL_width*FL_thickness;%[m3]
Ku_tmp=Ea/FL_volume;%[J/m3]
Ku=Ku_tmp*10;%[erg/cm3]<--[J/m3]
Hk_tmp2=2*Ku/Ms;%[Oe]
%Hk_tmp=Hk_tmp2*1e-4;%[tesla]<--[Oe]
Hk_tmp=500e-4;
Hk=[0,0,Hk_tmp];
Hext_tmp=800*1e-4;%tesla
Hext=[0,Hext_tmp,0];%tesla

runtime=10e-9;%[ns]
ts=1e-12;%[ps]
num_step=runtime/ts;
isposi=1;

theta_init=0.1/180*pi;

phi_init=0;
mx=sin(theta_init)*cos(phi_init);
my=sin(theta_init)*sin(phi_init);
mz=cos(theta_init);
m_init=[mx,my,mz];
p=[0,0,1];
Ie=0;

[mmx,mmy,mmz,tt,Icri]=rk4_4llg_PMA(ts,num_step,m_init,Ms,Hk,Hext,...
    alpha,P1,P2,p,Ie,FL_length*1e9,FL_width*1e9,FL_thickness*1e9);
toc

plot(tt,mmx,tt,mmy,tt,mmz,'linewidth',2);
legend('mx','my','mz')
ylim([-1,1])
xlabel('time(ns)','fontsize',15);ylabel('mxyz','fontsize',15)
set(gca,'fontsize',15)
print('zerostt', '-dpng', '-r300'); %<-Save as PNG with 300 DPI












