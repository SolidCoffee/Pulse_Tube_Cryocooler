%%%%%Pulse Tube Modeling Code%%%%%

%%%Initialization%%%
clc
clear

%%%Main Body Code%%%

%values
Tm = 300;           %Mean temperature, Tm, Kelvin
Dpt = 0.1;          %Diameter of Pulse Tube, Dpt, m
Lpt = .5;            %Length of Pulse Tube, Lpt,m
R = 2.0769;         %Heium Gas constant, R, 

P0 = 20;            %Ambient Pressure of the Pulse Tube, P0, Atmospheres
P0 = P0*101.325;    %Ambient Pressure of the Pulse Tube, P0, Kpa
P1 = 1;             %Oscillating Pressure from compressors, P1, Kpa

Th = 400;           %Temperature at the hot end of the pulse tube, Th, Kelvin
Tc = 150;           %Temperature at the cold end of the pulse tube, Tc, Kelvin

f = 45;             %Frequency of Oscillation, f, Hz
omega = 2*pi*f;     %Angular frequency, omega, Hz

%Initial calculations
Acpt = (pi/4)*Dpt^2;   %Cross sectional Area of Pulse Tube, Acpt, m^2
Vpt = Acpt*Lpt;        %Volume of Pulse Tube, Vpt, m^3


%Properties table A.4 from book
table = textread('helium.txt');

%Breaking out the table into a vector for each variable
T = table(:,1);
rho = table(:,2);
cp = table(:,3);
mu = table(:,4)*10^-7;
nu = table(:,5)*10^-6;
k = table(:,6)*10^-3;
alpha = table(:,7)*10^-6;
Pr = table(:,8);

%interpolating from table
Rhom = interp1(T,rho,Tm);


%symbolic math

syms P t mdoth mdotc

P = P0 + P1*cos(omega*t);

mdoth = (Th/Tc)*P1*cos(omega*t);
mdotc = mdoth + ((omega*Vpt)/(R*Tm))*P1*cos(omega*t + pi/2);

figure
fplot([mdoth, mdotc]);