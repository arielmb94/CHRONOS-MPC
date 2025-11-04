nx = 6; %[il wl ir wr]
nw = 5;
nu = 2;
nz = 6;
ny = 3;

A = [-Ra_l/La_l -K_l/La_l 0 0 0 0;
     K_l/J_l -F_l/J_l 0 0 0 0;
     0 0 -Ra_r/La_r -K_r/La_r 0 0;
     0 0 K_r/J_r -F_r/J_r 0 0;
     zeros(2,6)];    
B = [1/La_l 0;0 0;0 1/La_r;0 0;zeros(2,2)];
Bw = [zeros(1,nw);-0.01/J_l zeros(1,nw-1);zeros(1,nw);0 -0.01/J_r zeros(1,nw-2);zeros(2,nw)];
C = [0 1 0 0 1 0;
     0 0 0 1 0 1;
     1 0 1 0 0 0];
Cobs = C;
Dw = [zeros(ny,nw-3) diag([5000 5000 10])];%10 10 5
D = zeros(3,2);
sys = ss(A,[Bw B],C,[Dw D]);
ts = 0.02;
sysd = c2d(sys,ts,'Tustin');


%%

A = sysd.A;
Aobs = A;
Bw = sysd.B(:,1:nw);
Bu = sysd.B(:,nw+1:nw+nu);
Buobs = Bu;
Cz = diag([1 1 1 1 100 100]);
Dzw = zeros(nz,nw);

%%
A = [-Ra_l/La_l -K_l/La_l 0 0;
     K_l/J_l -F_l/J_l 0 0;
     0 0 -Ra_r/La_r -K_r/La_r;
     0 0 K_r/J_r -F_r/J_r];    
B = [1/La_l 0;0 0;0 1/La_r;0 0];
Bw = [0 0;-1/J_l 0;0 0;0 -1/J_r];
C = [0 1 0 0;
     0 0 0 1;
     1 0 1 0];
Gmotor = ss(A,[Bw B],C,0);
Gmotor.InputDelay = [0 0 tau_m tau_m];
