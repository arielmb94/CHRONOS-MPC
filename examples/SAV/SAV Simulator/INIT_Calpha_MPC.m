%%% Load controllers and sampling time
clear all
close all

load('xmotor.mat')
xmotor = x;

La_l = xmotor(1);
Ra_l = xmotor(2);
K_l = xmotor(3);
J_l = xmotor(4);
F_l = xmotor(5);

La_r = xmotor(1)
Ra_r = xmotor(2)
K_r = xmotor(3)
J_r = xmotor(4)
F_r = xmotor(5)

tau_m = xmotor(6)

load('xNL.mat')

C_sigma_2 = x(1)
C_sigma_1 = x(2)
C_sigma_0 = x(3)
C_alpha_f2 = x(4)
C_alpha_f1 = x(5)
C_alpha_f0 = x(6)
C_alpha_r2 = x(7)
C_alpha_r1 = x(8)
C_alpha_r0 = x(9)
eps = x(10)
wn = x(11)
tau = x(12)
Iz = x(13)
R = x(14)

servoSIM = ss([0 1; -wn^2 -2*eps*wn],[0; wn^2],[1 0],0);
servoSIM.InputDelay = tau;

Ts=1/50;

Observer_Motor_PI

Servo = ss([0 1; -wn^2 -2*eps*wn],[0; wn^2],[1 0],0);
d = 9; %tau/Ts
z = tf('z',Ts);
delay = ss(z^-d);

Servo = series(delay,c2d(Servo,Ts));

load("controllers.mat")

sav_lat_controller_init
%% Map

Cartoon

%% Inital conditions

global pathindex
pathindex = 1;
X0 = Track.X(1);
Y0 = Track.Y(1);
Psi0= -pi/2;
v0 = 1; %m/s
tp = 0.95;

vref = 1.0;


%% Run simulation and plot trajectory
pathindex = 1;


out = sim('CAR_SIM_MPC.slx',40);
close all
figure
plot(Track.X,Track.Y)
hold on
plot(out.X,out.Y)
axis equal

%%
figure
plot3(Track.X,Track.Y,ones(1,length(Track.X)))
hold on
plot3(Track.X,Track.Y,1.45*ones(1,length(Track.X)))
plot3(out.X,out.Y,out.Vx) 
axis equal
grid on
