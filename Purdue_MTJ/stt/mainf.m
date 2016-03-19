%% STT-_example.sp
clear all;clc
%% dimension
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
Hk_tmp=Hk_tmp2*1e-4;%[tesla]<--[Oe]
Hk=[0,0,Hk_tmp];
Hext_tmp=0*1e-4;%tesla
Hext=[0,0,Hext_tmp];%tesla

%% initial condition
isparallel=1;
switch isparallel
    case 0
        theta_init=179.9*pi/180;
mz=cos(theta_init); %0.1 degree tilt.
%m_init=[sqrt((1-mz^2)/2),sqrt((1-mz^2)/2),mz];
%m_init=[0,sqrt(1-mz^2),mz];
m_init=[sqrt(1-mz^2),0,mz];
Ie=100e-6;
    case 1
theta_init=0.1*pi/180;
mz=cos(theta_init); %0.1 degree tilt.
%m_init=[sqrt((1-mz^2)/2),sqrt((1-mz^2)/2),mz];
%m_init=[0,sqrt(1-mz^2),mz];
m_init=[sqrt(1-mz^2),0,mz];
Ie=-100e-6;
end

t1=3e-9;%s
num_step=3000;%number of steps
ts=t1/num_step;%time step size
p=[0,0,1];

[mmx,mmy,mmz,tt,Icri]=rk4_4llg_PMA(ts,num_step,m_init,Ms,Hk,Hext,...
    alpha,P1,P2,p,Ie,FL_length*1e9,FL_width*1e9,FL_thickness*1e9);
plot(tt,mmx,tt,mmy,tt,mmz,'linewidth',2);
legend('mx','my','mz')
xlabel('time(ns)','fontsize',15);ylabel('mxyz','fontsize',15)
set(gca,'fontsize',15)
%print('ptoap', '-dpng', '-r300'); %<-Save as PNG with 300 DPI


