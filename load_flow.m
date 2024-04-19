function [P12,P21,P13,P31,P24,P42,P34,P43, Losses] = load_flow(Linedata,Mag_V,Angle_V)
%Load flow & losses
B=Linedata(:,8)*sqrt(-1);%Multiplying the shunt element by j
y=Linedata(:,5)+(Linedata(:,6))*sqrt(-1);%Calculating the admittance (y) of each transmission line
V = Mag_V(:,1).*cosd(Angle_V(:,1))+sqrt(-1)*Mag_V(:,1).*sind(Angle_V(:,1));%Changing the voltage from phasor form to polar form
%==================================================Power Flow====================================================
P12 = 100*(conj(y(1)*(V(1)-V(2))+B(1)*V(1)));P21 = 100*(conj(y(1)*(V(2)-V(1))+B(1)*V(2))*V(2));
P13 = 100*(conj(y(2)*(V(1)-V(3))+B(2)*V(1))*V(1));P31 = 100*(conj(y(2)*(V(3)-V(1))+B(2)*V(3))*V(3));
P24 = 100*(conj(y(3)*(V(2)-V(4))+B(3)*V(2))*V(2));P42 = 100*(conj(y(3)*(V(4)-V(2))+B(3)*V(4))*V(4));
P34 = 100*(conj(y(4)*(V(3)-V(4))+B(4)*V(3))*V(3));P43 = 100*(conj(y(4)*(V(4)-V(3))+B(4)*V(4))*V(4));
%================================================Losses Calculations=============================================
Losses(1,1) = P12+P21;
Losses(2,1) = P13+P31;
Losses(3,1) = P24+P42;
Losses(4,1) = P34+P43;
end
