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
g = 9.81;
mu = 1;

%C_sigma;  %front cornering stiffnes
    
%%%% Step 1: Define Grid Points, Parameter Variations and Sampling Time

Problem.NumVaryingParameters = 1;
Vx = 0.5:0.01:2;   %m/s 
VxRate = 1;     %m/s^2
Ts = 1/50;      %s
Problem.Ts = Ts;

%%%% Step 2: For loop to define the generalazied plant system, basis
%%%% functions and Polytope due to parameter variation at each grid point

f1 = 0.25;
We = ss(tf([0.5 2*pi*f1],[1 2*pi*f1*0.01]));
We = c2d(We,Ts,'Tustin');
M = 30;
f2 = 10*M;
Wu = ss(tf([1 2*pi*f2/M], [0.001 2*pi*f2]));
Wu = c2d(Wu,Ts,'Tustin');

Wudt = ss(1/80);

Int = ss(tf(1,[1 0]));
Int = c2d(Int,Ts);

GridPoint = 0;
for VxPoint = 1:length(Vx)
    
    %Increase Grid Point Counter
    GridPoint = GridPoint+1;
    
    %Vx value at the grid point
    vx = Vx(VxPoint);
    
    %Build Generalized Plant at the grid point

    C_sigma = C_sigma_2*vx^2+C_sigma_1*vx+C_sigma_0;

    A = [-2*C_sigma/(m*vx)];
    B = [2*C_sigma*R/(m*vx)];
    G = ss(A,B,1,0);
    G = c2d(G,Ts,'zoh');
         
    Gp{GridPoint}=G;

    systemnames = 'G Int We Wu Wudt';
    inputvar = '[r;n;d;u]';
    outputvar = '[We;Wu;Wudt]';
    input_to_We = '[r-G-1.0*n]';
    input_to_Wu = '[Int]';
    input_to_Wudt = '[u]';
    input_to_G = '[Int]';
    input_to_Int = '[u]';
    cleanupsysic = 'yes';
    P = sysic;

    Problem.P{GridPoint} = P; % Generalized Plant at GridPoint
    
    %Basis Function
    %Problem.Basis{GridPoint} = [1,vx,1/vx];
    Problem.Basis{GridPoint} = [1,vx,1/vx];

    %Parameter Variation Politope Basis
    Problem.BasisVariation{GridPoint} = [1,(vx+Ts*VxRate),1/(vx+Ts*VxRate);
                                         1,(vx-Ts*VxRate),1/(vx-Ts*VxRate)];
    
end

Problem.TotalPoints = GridPoint;

%% Problem Resolution

[Klong,X,G,Y,CL,gamma] = lmiHinfExtSFDiscreteGrid2(Problem,1,15,'sdpt3')

%%

[Klong,X,Gk,CL,gamma] = lmiHinfExtSFDiscreteGrid_PL_V2(Problem,1,15,'sdpt3')
Welong = We;
Wulong = Wu;
pathindex = 1;


%%
close all
figure(1)
bodemag(1/We,'k')
figure(2)
bodemag(1/Wu,'k')
figure(3)
bodemag(1/Wudt,'k')

for i = 1:15:Problem.TotalPoints
    figure(1)
    hold on
    bodemag(CL{i}(1,1)/We)
    figure(2)
    hold on
    bodemag(CL{i}(2,1)/Wu)
    figure(3)
    hold on
    bodemag(CL{i}(3,1)/Wudt)
%     figure(3)
%     hold on
%     bodemag(CL{i}(2,1)/(Wu*Int))

end

