clc 
clear
buses = 4;%Number of buses
[Linedata,Ybus_matrix] = Ybus(buses);%Function that calculates the Ybus matrix
[Mag_V,Angle_V,iterations,PGen,QGen,PLoad,QLoad,Pslack,Qslack,Q4]=Calculate(buses,Ybus_matrix);
[P12,P21,P13,P31,P24,P42,P34,P43, Losses] = load_flow(Linedata,Mag_V,Angle_V);
 display_analysis(Mag_V,Angle_V,PGen,QGen,Pslack,Qslack,Q4,PLoad,QLoad,P12,P21,P13,P31,P24,P42,P34,P43,Losses,iterations);
 