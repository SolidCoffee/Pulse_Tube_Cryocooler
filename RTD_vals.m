%% RTD temp resistance table generator
clc
clear

R_0=1000;
a = 3.9083*10^-3;
b = -5.775*10^-7;
c = -4.183*10^-12;
i=1;


for T= -200:0
    R(i)=R_0*(1+a*T+b*T^2+c*(T-100)*T^3);
    i=i+1;
end
Table=[R']
writematrix(Table,'A.xls');

%% RTD temp resistance table generator
clc
clear


R_0=1000;
a = 3.9083*10^-3;
b = -5.775*10^-7;
c = -4.183*10^-12;
i=1;
for T= 0:850
    R(i)=R_0*(1+a*T+b*T^2);
    i=i+1;
end
Table=[R']
writematrix(Table,'B.xls');
