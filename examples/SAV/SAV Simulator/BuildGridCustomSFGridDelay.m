clear Problem

%%%% Step 0: System Parameters definition

%Car Parameters
%load('parametersCar.mat')

m=1.1934; % mass of the vehicle
l_r=0.1049; % distance from center of gravity to rear wheel
L=0.174;
l_f=L-l_r; 
Wi=0.087; % width of the car
tr = Wi/2;
%R = 0.065/2; %  Wheel radius
g = 9.81;
mu = 1;

%Parameters Identified
%Iz = 0.005;%0.005;
%C_f = 4.8438;%;2.1; %front cornering stiffnes
%C_r = 11.2442;%4.3; %rear cornering stiffnes

% Iz = 0.034;%0.005;
%C_f = C_alpha_f;%;2.1; %front cornering stiffnes
%C_r = C_alpha_r;%4.3; %rear cornering stiffnes

%%%% Step 1: Define Grid Points, Parameter Variations and Sampling Time

Problem.NumVaryingParameters = 1;

Vx = 0.6:0.3:1.8;   %m/s 
Problem.VaryingParam{1} = Vx;
VxRate = 1;     %m/s^2

Ts = 1/50;      %s
Problem.Ts = Ts;

%%%% Step 2: For loop to define the generalazied plant system, basis
%%%% functions and Polytope due to parameter variation at each grid point

Servo = ss([0 1; -wn^2 -2*eps*wn],[0; wn^2],[1 0],0);
d = 9; %tau/Ts
z = tf('z',Ts);
delay = ss(z^-d);

Servo = series(delay,c2d(Servo,Ts));


GridPoint = 0;
for VxPoint = 1:length(Vx)
    
    %Vx value at the grid point
    vx = Vx(VxPoint);
    %Increase Grid Point Counter
    GridPoint = GridPoint+1;
    Problem.Grid{GridPoint} = [vx];
    
    %Build Generalized Plant at the grid point

    C_f = C_alpha_f2*vx^2+C_alpha_f1*vx+C_alpha_f0;
    C_r = C_alpha_r2*vx^2+C_alpha_r1*vx+C_alpha_r0;
    
    A = [-(C_f+C_r)/(m*vx) -vx-(C_f*l_f-C_r*l_r)/(m*vx);
     -(C_f*l_f-C_r*l_r)/(Iz*vx) -((C_f*l_f^2)+(C_r*l_r^2))/(Iz*vx)];
    B = [C_f/m; C_f*l_f/Iz];
    C = [0 1];
    G = ss(A,B,C,0);
    G = c2d(G,Ts,'zoh');

    f2 = 1;%10                                                      %%5
    Wu = ss(tf([1 2*pi*f2/0.25], [0.0001 2*pi*f2])); %0.5 0.0001      %%0,4
    Wu = c2d(Wu,Ts,'Tustin');

    M = 2;
    f1 = 0.5;
    Ae = 0.001;
    We = ss(tf([1/M 2*pi*f1],[1 2*pi*f1*Ae]));
    We = c2d(We,Ts,'Tustin');

    dist = ss(vx^2);

    systemnames = 'G Servo We Wu dist';
    inputvar = '[desR(1); d; n;delta(1)]';
    outputvar = '[We; Wu]';
    input_to_Wu = '[delta]';
    input_to_We = '[desR-G-1*n]';
    input_to_G = '[Servo+0.1*dist]';
    input_to_Servo = '[delta]';
    input_to_dist = '[d]';
    cleanupsysic = 'yes';
    P = sysic;
    Problem.P{GridPoint} = P; % Generalized Plant at GridPoint

    %Basis Function
    %Problem.Basis{GridPoint} = [1,vx,1/vx];
    Problem.Basis{GridPoint} = [1,vx,vx^2,1/vx];

    %Parameter Variation Politope Basis
    Problem.BasisVariation{GridPoint} = [1,vx+Ts*VxRate,(vx+Ts*VxRate)^2,1/(vx+Ts*VxRate);
        1,vx-Ts*VxRate,(vx-Ts*VxRate)^2,1/(vx-Ts*VxRate)];
    
end

Problem.TotalPoints = GridPoint;

%% Problem Resolution
%[Kval,Xval,Gval,CL,gopt] = lmiHinfExtSFDiscreteGrid_PL_V2(Problem,1,15,'sdpt3')

%%

[Klat,X,G,Y,CL,gopt] = lmiHinfExtSFDiscreteGrid2(Problem,1,15,'sdpt3')
Welat = We;
Wulat = Wu;

pathindex = 1;

%%
close all
figure(1)
bodemag(1/We,'k')
figure(2)
bodemag(1/Wu,'k')

for i = 1:Problem.TotalPoints
    figure(1)
    hold on
    bodemag(CL{i}(1,1)/We)
    figure(2)
    hold on
    bodemag(CL{i}(2,1)/Wu)
%     figure(3)
%     hold on
%     bodemag(CL{i}(1,2)/(We*Vx(i)^2))

end

