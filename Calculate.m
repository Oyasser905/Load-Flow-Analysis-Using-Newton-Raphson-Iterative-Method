function [Mag_V,Angle_V,iterations,PGen,QGen,PLoad,QLoad,Pslack,Qslack,Q4]=Calculate(buses,Ybus_matrix)
%1. Slack bus | 2. PQ bus | 3. PV bus
%           Bus | Pgen | Qgen | Pload | Qload |  V  | Delta  
%           no. |  MW  | MVAR |   MW  | MVAR  |  pu | Degree      
Busdata = [  1     0      0      50     30.99   1.00    0;
             2     0      0      170    105.35  1.00    0;
             3     0      0      200    123.94  1.00    0;
             4    318     0      80     49.58   1.02    0 ];    
         
Mag_V=Busdata(:,6);%Bus voltage magnitude
Angle_V=Busdata(:,7);%Bus voltage angle in degrees
PGen=Busdata(:,2);%Bus generated active power
QGen=Busdata(:,3);%Bus generated reactive power
PLoad=Busdata(:,4);%Bus load active power
QLoad=Busdata(:,5);%Bus load reactive power
S=100; %Base apparent power
tol_1 = 1;%Change of the magnitude of the voltage(difference between correct value and the current value)
tol_2 = 1;%Change of the angle of the voltage(difference between correct value and the current value)
no_pq=2;%Number of PQ buses
no_pv=1;%Number of PV buses
iterations = 0;
Mag_Ybus=abs(Ybus_matrix); %Magnitude of each element of Ybus
Angle_Ybus=rad2deg(angle(Ybus_matrix)); %Angle in degree of each element of Ybus in Degrees
%==================================================Iterations====================================================
while ((abs(tol_1)>10e-6) || (abs(tol_2)>10e-6))
%Calculated active power
Psum = 0;
for i = 2:buses
    for n = 1:buses
        if n~=i
            Pi=(Mag_V(i)*Mag_V(n)*Mag_Ybus(i,n))*cosd(Angle_Ybus(i,n)+Angle_V(n)-Angle_V(i));
        end
        Psum = Psum + Pi;
        Pi =0;
    end
    Pcalc(i-1)=((Mag_V(i))^2)*real(Ybus_matrix(i,i))+Psum;
    Psum = 0;
end
%Schedueled active power
for j = 2:buses
    Psch(j-1)=(PGen(j)-PLoad(j))/(S);
end
%Delta active power
deltaP=Psch-Pcalc;

%Calculated reactive power
Qsum = 0;
for i=2:buses-1
    for n=1:buses
        if n~=i
             Qi=(Mag_V(i)*Mag_V(n)*Mag_Ybus(i,n))*sind(Angle_Ybus(i,n)+Angle_V(n)-Angle_V(i)); 
        end
        Qsum = Qsum + Qi;
        Qi =0;
    end
    Qcalc(i-1)=-1*((Mag_V(i))^2)*imag(Ybus_matrix(i,i))-Qsum;
    Qsum = 0;
end

%Schedueled reactive power
for j = 2:buses-1
    Qsch(j-1)=(QGen(j)-QLoad(j))/(S);
end

%Delta reactive power
deltaQ=Qsch-Qcalc;

%Mismatch matrix
mis_m=zeros(5,1);
for i=1:buses-1
    mis_m(i)=deltaP(i);
end
for i=1:no_pq
    mis_m(i+buses-1)=deltaQ(i);
end
%==================================================Jacobian Matrix Calculations====================================================
J_dim=2*buses-2-no_pv;
J=zeros(J_dim,J_dim);
%=================================================================J11==============================================================
J11=zeros(buses-1,buses-1);
for i = 2:buses                      
    for j=2:buses
        if i~=j
            J11(i-1,j-1)=-(Mag_V(i)*Mag_V(j)*Mag_Ybus(i,j))*sind(Angle_Ybus(i,j)+Angle_V(j)-Angle_V(i));
        elseif (i==j && i~=4)
            J11(i-1,j-1)=-Qcalc(i-1)-((Mag_V(i))^2)*imag(Ybus_matrix(i,i));
        elseif(i==j && i==4)
             J11(i-1,i-1)=-J11(i-1,i-2)-J11(i-1,i-3);
        end
        J(i-1,j-1)=J11(i-1,j-1);
    end
end
%================================================================J21===================================================================
J21=zeros(no_pq,no_pq+no_pv);
for i = 2:buses-1                    
    for j=2:buses
        if i~=j
            J21(i-1,j-1)=-(Mag_V(i)*Mag_V(j)*Mag_Ybus(i,j))*cosd(Angle_Ybus(i,j)+Angle_V(j)-Angle_V(i));
        elseif (i==j && i~=4)
            J21(i-1,j-1)=Pcalc(i-1)-((Mag_V(i))^2)*real(Ybus_matrix(i,i));
        end
        J(i+2,j-1)=J21(i-1,j-1);
    end
end
%==============================================================J12====================================================================
J12=zeros(no_pq,no_pq);
for i = 2:buses                  
    for j=2:buses-1  
        if i~=j
            J12(i-1,j-1)=(Mag_V(i)*Mag_V(j)*Mag_Ybus(i,j))*cosd(Angle_Ybus(i,j)+Angle_V(j)-Angle_V(i));
        elseif (i==j)
            J12(i-1,j-1)=Pcalc(i-1)+((Mag_V(i))^2)*real(Ybus_matrix(i,i));
        end
        J(i-1,j+2)=J12(i-1,j-1);
    end
end
%============================================================J22======================================================================
J22=zeros(no_pq,no_pq);
for i = 2:buses-1                      
    for j=2:buses-1
        if i~=j
            J22(i-1,j-1)=-(Mag_V(i)*Mag_V(j)*Mag_Ybus(i,j))*sind(Angle_Ybus(i,j)+Angle_V(j)-Angle_V(i));
        elseif (i==j && i~=4)
            J22(i-1,j-1)=Qcalc(i-1)-((Mag_V(i))^2)*imag(Ybus_matrix(i,i));
        end
        J(i+2,j+2)=J22(i-1,j-1);
    end
end
%==============Creating the updates matrix, updating the voltage (Magnitude & angle), & calculating the new tolerance==================
updates=inv(J)*mis_m;
for x=1:3 % Delta elements to calculated in degrees
    updates(x)=rad2deg(updates(x));
end
deltaV(3:4)=updates(3:4).*Mag_V(2:3);
deltaS = updates(1:3);
%Updating the voltage angles
Angle_V(2:4,:) = updates(1:3)+Angle_V(2:4);
Mag_V(2:3) = (updates(4:5).* Mag_V(2:3))+ Mag_V(2:3);
tol_1=max(abs(deltaV));
tol_2=max(abs(deltaS));
%=======================================================Slack Bus Powers==============================================================================
Psum = 0;
Qsum = 0;
for i = 2:buses
    P1=(Mag_V(1)*Mag_V(i)*Mag_Ybus(1,i))*cosd(Angle_Ybus(1,i)+Angle_V(i)-Angle_V(1));
    Q1=(Mag_V(1)*Mag_V(i)*Mag_Ybus(1,i))*sind(Angle_Ybus(1,i)+Angle_V(i)-Angle_V(1)); 
    Psum = Psum + P1;
    Qsum = Qsum + Q1;
end
Pslack=((Mag_V(1))^2)*real(Ybus_matrix(1,1))+Psum;
Qslack=-1*((Mag_V(1))^2)*imag(Ybus_matrix(1,1))-Qsum;   
%Q for bus 4 (PV bus)
Qsum = 0;
for i = 1:buses-1
    Q4=(Mag_V(4)*Mag_V(i)*Mag_Ybus(4,i))*sind(Angle_Ybus(4,i)+Angle_V(i)-Angle_V(4)); 
    Qsum = Qsum + Q4;
end
Q4=-1*((Mag_V(4))^2)*imag(Ybus_matrix(4,4))-Qsum;
%Counting the number of iterations
iterations = iterations + 1;
end
Pslack = (Pslack * S) + PLoad(1); 
Qslack = (Qslack * S) + QLoad(1);
Q4 = (Q4 * S) + QLoad(4);
end 
