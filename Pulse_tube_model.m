%%%%%Pulse Tube Modeling Code - Active Piston displacer%%%%%




%%%%%REFERENCES%%%%%

%%%1%%%
%MAIN REFERENCE
%https://cryocooler.org/resources/Documents/C19/201.pdf#page=1

%ADDITIONAL REFERENCE (ONLY USED AS GUIDANCE AND SANITY CHECK OF MAIN
%REFERENCE)
%https://trc.nist.gov/cryogenics/Papers/Pulse_Tube_Cryocoolers/2010-Development_of_Miniature_High_Frequency_PTC.pdf

%%%Heat Transfer TextBook%%%
%FOR PROPERTIES TABLE
%Fundamentals_of_Heat_and_Mass_Transfer_7th_Edition_Incropera



%%%Initialization%%%
clc                 %Clears Command Window
clear               %Clears Workspace Window (Variables in Memory)



%%%User defined values%%%

%Geometry
Dpt = 0.010;                %Diameter of Pulse Tube, Dpt, m
Lpt = .1;                   %Length of Pulse Tube, Lpt,m

%Other given dimensions
Vregen = 20*10^-6;          %Volume of Regenerator, Vregen, m^3*10^-6 OR cc
Vcd = 20*10^-6;             %Volume of dead volume in the compressor, Vcd, m^3*10^-6 OR cc
Vco = Vcd;                  %Volume of the compressor, Vco, m^3*10-6 or cc, keep as Vcd for now
Vap = Vco;                  %Swept Volume of the active piston, Vco, m^3*10^-6 or cc, keep as Vco for now.
Vhx = 2*10^-6;              %Volume of the heat exchangers, Vhx, m^3*10^-6 OR cc
Vac = Vhx; Vhhx = Vhx; Vchx = Vhx;     %Heat exchangers have same volume (Both cold and Hot and after cooler)

%Pressures
P0 = 20;            %Ambient Pressure of the Pulse Tube, P0, Atmospheres
                    %Static Pressure of the fluid within the tube when the
                    %compressor is off, read by pressure transducer

P1 = 1;             %Oscillating Pressure from compressors, P1, Atmospheres
                    %Value of the oscillating pressure wave of the fluid caused by the compressor
                    %Value is a peak value, aka 1 atm peak (Not peak to
                    %peak) value would be read by pressure transducers

%Temperatures
Th = 300;           %Temperature at the hot end of the pulse tube, Th, Kelvin
Tc = 100;           %Temperature at the cold end of the pulse tube, Tc, Kelvin
Tac = 310;          %Temperature at the after cooler heat exchanger, Tac, Kelvin
Tco = 320;          %Temperature at the compressor, Tco, Kelvin
                    %Values would be measured using thermocouples positioned at each end

%Frequency
f = 45;             %Drive frequency, Frequency of Oscillation, f, Hz



%Constants and properties
C1 = 6*10^-10;     %Flow Coefficient, Get from page 8 (page 306) reference document
gamma = 1.67;      %Specific heat ratio of helium, gamma, dimensionless (cp/cv)
R = 2076.9;        %Heium Gas constant, R, J/Kg*K or pa*m^3 / Kg*K



%%%Main Body Code%%%


%Intermediate Calculations
P0 = P0*101.325*1000;    %Ambient Pressure of the Pulse Tube, P0, pa
                         %Converts from atmosphere to pa
                        
P1 = P1*101.325*1000;    %Oscillating Pressure from compressors, P1, pa
                         %Converts from atmosphere to pa

Tm = (Th + Tc)/2;                    %Arithmatic Mean temperature, Tm, Kelvin
                                     %Would be a mean value of the temp measured from the
                                     %thermocouple in the flow of the fluid

Tlm = (Th - Tc)/(log(Th/Tc));        %Log mean temperature, Tlm, Kelvin

Acpt = (pi/4)*Dpt^2;   %Cross sectional Area of Pulse Tube, Acpt, m^2

Vpt = Acpt*Lpt;        %Volume of Pulse Tube, Vpt, m^3

omega = 2*pi*f;     %Angular frequency, omega, rad/s


%Properties of Helium at different temperatures
%table A.4 from Heat transfer book, page 997
%Fundamentals_of_Heat_and_Mass_Transfer_7th_Edition_Incropera
table = textread('helium.txt');

%Breaking out the table into a vector for each variable
T = table(:,1);                     %Temperature, T, Kelvin
rho = table(:,2);                   %Density, rho, Kg/m^3
cp = table(:,3);                    %Specific heat at constant pressure, Cp, KJ/Kg*K
mu = table(:,4)*10^-7;              %Dynamic Viscosity, Mu, N*s/m^2
nu = table(:,5)*10^-6;              %Kinematic Viscosity, Nu, m^2/s
k = table(:,6)*10^-3;               %Thermal conductivity, k, W/m*K
alpha = table(:,7)*10^-6;           %Thermal Diffusivity, alpha, m^2/s
Pr = table(:,8);                    %Prandtl Number, Pr, Dimensionless

%interpolating from table
Rhom = interp1(T,rho,Tm);           %Mean value of density (rho), rhom, Kg/m^3
    %Value is considered a mean value because it is interpolated from the
    %table at the mean temperature value Tm




%%%Mass flow calc%%%

%all to be treated as phasors
%consider as stationary phasors, aka vector



%Going from bottom to top
%Magnitude of the mass flow vector at the piston
map_mag = (P0*omega*Vap)/(4*R*Tc);



%Vertical component of the mass flow vector at mc
mc_vert = -((omega*P1*Vregen)/(R*Tlm))/2;

%Vertical component of the mass flow vector at mchx
mchx_vert = -((omega*P1*Vchx)/(R*Tc)) + mc_vert;

%Vertical component of the mass flow vector mh
mh_vert = -((omega*P1*Vpt)/(gamma*R*Tc)) + mchx_vert;

%Vertical component of the mass flow vector mp
mp_vert = -((omega*P1*Vhhx)/(R*Tc)) + mh_vert;

%Vertical component of the mass flow vector at the active piston
map_vert = -(omega*P1*Vap)/(2*R*Tc) + mp_vert;



%Phase angle of the active piston mass flow relative to the middle of the
%regenerator
phi_ap = asind(map_vert / map_mag);

%Horizontal component of the active piston mass flow vector, which is also the
%mass flow at the middle of the regenerator. This is the minimum mass flow
%in the system
map_hor = map_mag*cosd(phi_ap);



%Vertical component of the mass flow vector at the aftercooler
mac_vert = ((omega*P1*Vregen)/(R*Tlm))/2;

%Vertical component of the mass flow vector at m1
m1_vert = ((omega*P1*Vac)/(R*Tac)) + mac_vert;

%Vertical component of the mass flow vector at mco
mco_vert = (P1*omega*Vco)/(2*R*Tco) + m1_vert;












%Vector of mass flow at compressor
mco_vec = [map_hor, mco_vert];
mco_mag = norm(mco_vec)

%Vector of m1 mass flow vector
m1_vec = [map_hor, m1_vert];
m1_mag = norm(m1_vec);

%Vector of mass flow at aftercooler
mac_vec = [map_hor, mac_vert];
mac_mag = norm(mac_vec);

%Vector of mass flow at middle of regenerator
mrm_vert = 0;
mrm_vec = [map_hor, mrm_vert];
mrm_mag = norm(mrm_vec);

%Vector of mass flow at mc
mc_vec = [map_hor, mc_vert];
mc_mag = norm(mc_vec);

%Vector of mass flow at chx
mchx_vec = [map_hor, mchx_vert];
mchx_mag = norm(mchx_vec);

%Vector of mass flow mh
mh_vec = [map_hor, mh_vert];
mh_mag = norm(mh_vec);

%Vector of mass flow mp
mp_vec = [map_hor, mp_vert];
mp_mag = norm(mp_vec);

%Vector of mass flow at active piston
map_vec = [map_hor, map_vert];
map_mag = norm(map_vec);



%Phase angles
phi_1 = atand(m1_vert / map_hor);
phi_ac = atand(mac_vert / map_hor);
phi_rm = atand(mrm_vert / map_hor);
phi_c = atand(mc_vert / map_hor);
phi_chx = atand(mchx_vert / map_hor);
phi_h = atand(mh_vert / map_hor);
phi_p = atand(mp_vert / map_hor);



%Phase angle of the compressor mass flow relative to the middle of the
%regenerator
phi_co = atand(mco_vert / map_hor);




%In radians
theta_co = atan(mco_vert / map_hor);
theta_1 = atan(m1_vert / map_hor);
theta_ac = atan(mac_vert / map_hor);

theta_c = atan(mc_vert / map_hor);
theta_chx = atan(mchx_vert / map_hor);
theta_h = atan(mh_vert / map_hor);
theta_p = atan(mp_vert / map_hor);
theta_ap = atan(map_vert / map_hor);


%Angle between compressor and active piston
alpha = phi_co - phi_ap;

%%%%%symbolic math%%%%%

%all to be treated as sinusoidal waves
    %This is an equivalent represenation of a cyclic phasor (Not a vector)
        %A vector is a stationary phasor

syms t        %Creating symbolic variables
    %Syms are treated as different objects, as syms, and therefore need to
    %be handled and considered as such
   
P = P0 + P1*cos(omega*t);
%Pressure waveform, P, Kpa
%Equation 1 from reference 1
%Symbolic waveform of the pressure wave in the system
    %Static component P0, dynamic component from the compressor P1 driven
    %at compressor frequency

    
%Top Starting point
mco = mco_mag*cos(omega*t + theta_co);

%M1
m1 = m1_mag*cos(omega*t + theta_1);

%Mac
mac = mac_mag*cos(omega*t + theta_ac);

%Mrm
mrm = mrm_mag*cos(omega*t);

%Mc
mc = mc_mag*cos(omega*t + theta_c);

%Mchx
mchx = mchx_mag*cos(omega*t + theta_chx);

%Mh
mh = (mh_mag / (Th/Tc))*cos(omega*t + theta_h);

%Mp
mp = (mp_mag / (Th/Tc))*cos(omega*t + theta_p);

%Map
map = map_mag*cos(omega*t + theta_ap);




    

%Plotting of the two mass flow rates
    %You will need to zoom in to see the phase shift, its small but its
    %there

%Plot of all mass flows vs time for one cycle 
figure
graph = fplot([mco, m1, mac, mrm, mc, mchx, mh, mp, map], [0 0.022]);
graph(1).Color = 'r';
graph(2).Color = 'k';
graph(3).Color = 'g';
graph(4).Color = 'c';
graph(5).Color = 'm';
title('Mass flow rate Vs. Time, waveforms, one cycle')
xlabel('Time, t, seconds')
ylabel('Mass flow rate, mdot, Kg/s')
legend('Mdot co: compressor top','Mdot 1: Compressor bottom, AC top', 'Mdot ac: AC bot, RG top', 'Mdot rm: Regen middle, Mdot c: RG bot, CHX top', 'Mdot CHX: CHX bot, PT top', 'RG top, AC bot', 'Mdot h top: PT bot, hhx top', 'Mdot p: hhx bot, active piston top', 'Mdot ap: active piston bottom')

%Displaying the phase angles in the command window


phi_co
phi_1
phi_ac
phi_rm
phi_c
phi_chx
phi_h
phi_p
phi_ap

alpha



% %Plot of mass flow across CHX vs mass flow at orifice, for only one cycle
% figure
% fplot([mo, mc], [0 0.022])
% title('Mass flow rate Vs. Time, waveforms, one cycle')
% xlabel('Time, t, seconds')
% ylabel('Mass flow rate, mdot, Kg/s')
% legend('Mdot orifice end', 'Mdot chx end')



%Cooling power Qdotc
theta = abs(abs(theta_ap) - abs(theta_chx));

Qdotc = ((R*Tc*P1*mc)/(2*P0))*cos(theta);

figure
fplot([Qdotc], [0 0.04])
title('Cooling power Vs. Time, waveform, one cycle')
xlabel('Time, t, seconds')
ylabel('Cooling power, Qdot, W')
legend('Cooling power at CHX')



%Table combining all calculated values of the phasor analysis
A = [map_hor, map_hor, map_hor, map_hor, map_hor, map_hor, map_hor, map_hor, map_hor,; 
    mco_vert, m1_vert, mac_vert, mrm_vert, mc_vert, mchx_vert, mh_vert, mp_vert, map_vert;
    mco_mag, m1_mag, mac_mag, mrm_mag, mc_mag, mchx_mag, mh_mag, mp_mag, map_mag;
    phi_co, phi_1, phi_ac, phi_rm, phi_c, phi_chx, phi_h, phi_p, phi_ap];


%Exporting the table to excel so that data can be analyzed and graphed
%within excel
xlswrite('phasor_analysis.xlsx',A)



%Function waveforms
    %By removing the semicolon, this will display the full functional
    %description of the waveform
        %AKA, it will tell you what the exact sine function is, which can
        %be plugged into wolfram alpha or elsewhere

%Mass flow rate at CHX
mc
%Cooling power at CHX
Qdotc
