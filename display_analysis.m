function display_analysis(Mag_V,Angle_V,PGen,QGen,Pslack,Qslack,Q4,PLoad,QLoad,P12,P21,P13,P31,P24,P42,P34,P43,Losses,iterations)
%Display System Analysis
fprintf('Number of iterations = %d\n',iterations)
fprintf('X---------------------------------Bus information--------------------------------X----------------------------------Line flow----------------------------------X\n')
fprintf('Bus \t  \t Volts \t Angle \t\tX----Generation----X \t\t X----Load----X      Bus \t \tTo \t \t\t\t\t     Lineflow \t\t\t\t    Line Losses\n')
fprintf('no.   Name \t (p.u.)  (deg.)  \t(MW) \t\t\t(MVAR) \t \t(MW) \t  (MVAR)     type \t Bus Name \t\t\t\t (MW) \t\t (MVAR) \t\t (MW) \t\t (MVAR)\n')
fprintf('X--------------------------------------------------------------------------------------------------------------------------------------------------------------X\n')
fprintf('1 \t  Birch   %.2f \t %.3f \t\t %.3f \t    %.3f\t \t %d \t  %.2f \t SL \t\t2 \t\t  Elm\t\t %.3f \t %.3f \t\t %.3f \t\t %.3f\n',Mag_V(1),Angle_V(1),Pslack,Qslack,PLoad(1),QLoad(1),real(P12),imag(P12),real(Losses(1)),imag(Losses(1)))
fprintf('                                                                                            3 \t\t  Pine\t\t %.3f \t %.3f \t\t %.3f \t\t %.3f\n',real(P13),imag(P13),real(Losses(2)),imag(Losses(2)))
fprintf('2 \t  Elm     %.2f \t%.3f \t\t %d \t\t\t\t   %d \t \t %d \t  %.2f \t PQ \t\t1 \t\t  Birch\t\t %.3f \t %.3f \t\t %.3f \t\t %.3f\n',Mag_V(2),Angle_V(2),PGen(2),QGen(2),PLoad(2),QLoad(2),real(P21),imag(P21),real(Losses(1)),imag(Losses(1)))
fprintf('                                                                                            4 \t\t  Maple\t\t %.3f \t %.3f \t\t %.3f \t\t %.3f\n',real(P24),imag(P24),real(Losses(3)),imag(Losses(3)))
fprintf('3 \t  Pine    %.2f \t%.3f \t\t %d \t\t\t\t   %d \t \t %d \t  %.2f\t PQ \t\t1 \t\t  Birch\t\t %.3f \t %.3f \t\t %.3f \t\t %.3f\n',Mag_V(3),Angle_V(3),PGen(3),QGen(3),PLoad(3),QLoad(3),real(P31),imag(P31),real(Losses(2)),imag(Losses(2)))
fprintf('                                                                                            4 \t\t  Maple\t\t %.3f \t %.3f \t\t %.3f \t\t %.3f\n',real(P34),imag(P34),real(Losses(4)),imag(Losses(4)))
fprintf('4 \t  Maple   %.2f \t %.3f \t\t%d \t\t\t%.3f \t %d \t  %.2f \t PV \t\t2 \t\t  Elm\t\t %.3f \t %.3f \t\t %.3f \t\t %.3f\n',Mag_V(4),Angle_V(4),PGen(4),Q4,PLoad(4),QLoad(4),real(P42),imag(P42),real(Losses(3)),imag(Losses(3)))
fprintf('                                                                                            3 \t\t  Pine\t\t %.3f \t %.3f \t\t %.3f \t\t %.3f\n',real(P43),imag(P43),real(Losses(4)),imag(Losses(4)))
fprintf('X--------------------------------------------------------------------------------------------------------------------------------------------------------------X\n')
end