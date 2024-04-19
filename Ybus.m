function [Linedata,matrix] = Ybus(buses)
%         |  From |  To   |    R       |    X    |     G      |     B      |  Charging |  Y/2
%         |  Bus  | Bus   |    pu      |    pu   |     pu     |     pu     |    MVAR   |   pu 
Linedata=[    1      2       0.01008     0.05040   3.815629    -19.078144      10.25     0.05125;
              1      3       0.00744     0.03720   5.169561    -25.847809       7.75     0.03875;
              2      4       0.00744     0.03720   5.169561    -25.847809       7.75     0.03875;
              3      4       0.01272     0.06360   3.023705    -15.118528      12.75     0.06375];
Ybus_matrix=zeros(buses);               % Creation of 4x4 Ybus full of zeros
%Calculating diagonal elements
bus_1 = (Linedata(1,5) + Linedata(2,5)) + (Linedata(1,6) + Linedata(2,6) + Linedata(1,8) + Linedata(2,8))*i;
bus_2 = (Linedata(1,5) + Linedata(3,5)) + (Linedata(1,6) + Linedata(3,6) + Linedata(1,8) + Linedata(3,8))*i;
bus_3 = (Linedata(2,5) + Linedata(4,5)) + (Linedata(2,6) + Linedata(4,6) + Linedata(2,8) + Linedata(4,8))*i;
bus_4 = (Linedata(3,5) + Linedata(4,5)) + (Linedata(3,6) + Linedata(4,6) + Linedata(3,8) + Linedata(4,8))*i;
bus_i = [bus_1 bus_2 bus_3 bus_4];
for r = 1: buses
    for c = 1:buses
        if r == c
            Ybus_matrix(r,c) = bus_i(r);
        end
    end
end
%Calculating off-diagonal elements
Ybus_matrix(1,2) = (Linedata(1,5) + (Linedata(1,6)*i))*-1;
Ybus_matrix(1,3) = (Linedata(2,5) + (Linedata(2,6)*i))*-1;
Ybus_matrix(2,4) = (Linedata(3,5) + (Linedata(3,6)*i))*-1;
Ybus_matrix(3,4) = (Linedata(4,5) + (Linedata(4,6)*i))*-1;
%Skew symm. elements
Ybus_matrix(2,1) = Ybus_matrix(1,2);
Ybus_matrix(3,1) = Ybus_matrix(1,3);
Ybus_matrix(4,2) = Ybus_matrix(2,4);
Ybus_matrix(4,3) = Ybus_matrix(3,4);

matrix = Ybus_matrix;
end 