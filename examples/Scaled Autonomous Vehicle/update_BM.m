function A = update_BM(vx)

m=1.1934; % mass of the vehicle
l_r=0.1049; % distance from center of gravity to rear wheel
L=0.174;
l_f=L-l_r; 

Iz = 0.005717846957287;

C_alpha_f2 = -0.494189184123490;
C_alpha_f1 = 6.312178044425870;
C_alpha_f0 = -1.993719141364098;
C_alpha_r2 = 2.737538481248687;
C_alpha_r1 = 9.211804162490232;
C_alpha_r0 = -3.192821253237006;


C_f = C_alpha_f2*vx^2+C_alpha_f1*vx+C_alpha_f0;
C_r = C_alpha_r2*vx^2+C_alpha_r1*vx+C_alpha_r0;

A = [-(C_f+C_r)/(m*vx) -vx-(C_f*l_f-C_r*l_r)/(m*vx);
    -(C_f*l_f-C_r*l_r)/(Iz*vx) -((C_f*l_f^2)+(C_r*l_r^2))/(Iz*vx)];

end