%% DOF motor controller
Ts=1/50;

A = [-Ra_l/La_l -K_l/La_l; 
     K_l/J_l -F_l/J_l];

B = [1/La_l ; 0];

C = [0 1]; 

G = ss(A,B,C,0);
G = c2d(G,Ts,'Tustin');

n = round(tau_m/Ts);
z = tf('z',Ts);
Delay = ss(z^-(n));

f1 = 1;%0.5
We = ss(tf([0.5 2*pi*f1],[1 2*pi*f1*0.001]));%0.0001
We = c2d(We,Ts,'Tustin');
WeMotor = We;

Mu = 0.1;
f2 = 0.1;%0.1
Wu = ss(tf([1 2*pi*f2/Mu], [0.00001 2*pi*f2]));%0.01 // 0.0001 
Wu = c2d(Wu,Ts,'Tustin');
WuMotor = Wu;


%% DOF Hinfsyn

systemnames = 'G Delay We Wu';
inputvar = '[wref(1); n(1); V(1)]';
outputvar = '[We; Wu; wref-G-10*n]';
input_to_Wu = '[V]';
input_to_We = '[wref-G-10*n]';
input_to_G = '[Delay]';
input_to_Delay = '[V]';
cleanupsysic = 'yes';
P = sysic;

[Kbldc,CL,gamma] = hinfsyn(P,1,1,'method','ric')

%[Kbldc,gopt,gval,Pcl,CL,Xval,Yval] = lmiHinfDOFDiscrete_PL(P,1,1,0,15,'sdpt3')

close all
figure
bodemag(1/We,CL(1,1)/We)
figure
bodemag(1/Wu,CL(2,1)/Wu)
