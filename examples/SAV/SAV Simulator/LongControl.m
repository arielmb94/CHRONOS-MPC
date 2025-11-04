%Car Parameters
m=1.19; % mass of the vehicle
Ts = 1/50;

longModel = ss(0,1/m,1,0);
G = c2d(longModel,Ts);

f1 = 0.5;
We = ss(tf([0.5 2*pi*f1],[1 2*pi*f1*0.001]));
We = c2d(We,Ts,'Tustin');
f2 = 10;
Wu = ss(tf([1 2*pi*f2/10], [0.0001 2*pi*f2]));
Wu = c2d(Wu,Ts,'Tustin');

systemnames = 'G We Wu';
inputvar = '[r;u]';
outputvar = '[We;Wu;r-G]';
input_to_We = '[r-G]';
input_to_Wu = '[u]';
input_to_G = '[u]';
cleanupsysic = 'yes';
P = sysic



%%

[Klong,CL,gamma,info] = hinfsyn(P,1,1,'method','lmi')
figure
bodemag(1/We,CL(1,1)/We)
figure
bodemag(1/Wu,CL(2,1)/Wu)

% %%
% [Klong,gopt,gval,Xval,CL,R,S] = lmiHinfDOFDiscrete_PL(P,1,1,1,20,'mosek');
% figure
% bodemag(1/We,CL(1,1)/We)
% figure
% bodemag(1/Wu,CL(2,1)/Wu)