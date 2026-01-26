%%% Load controllers and sampling time

clear all
close all

% Sampling Time
Ts=1/50;

% Motor Parameters
load('xmotor.mat')
xmotor = x;

La_l = xmotor(1);
Ra_l = xmotor(2);
K_l = xmotor(3);
J_l = xmotor(4);
F_l = xmotor(5);

La_r = xmotor(1);
Ra_r = xmotor(2);
K_r = xmotor(3);
J_r = xmotor(4);
F_r = xmotor(5);

tau_m = xmotor(6);

% SAV parameters
load('xNL.mat')

C_sigma_2 = x(1);
C_sigma_1 = x(2);
C_sigma_0 = x(3);
C_alpha_f2 = x(4);
C_alpha_f1 = x(5);
C_alpha_f0 = x(6);
C_alpha_r2 = x(7);
C_alpha_r1 = x(8);
C_alpha_r0 = x(9);
eps = x(10);
wn = x(11);
tau = x(12);
Iz = x(13);
R = x(14);

% Steering Actuator Model
servoSIM = ss([0 1; -wn^2 -2*eps*wn],[0; wn^2],[1 0],0);
servoSIM.InputDelay = tau;

% Disturbance Observer for the Motor Control
Init_Observer_Motor_PI

% Longitudinal and Motor Controllers
load("controllers.mat")

% Init MPC for the lateral controller
sav_lat_controller_init
%% Load Map

Cartoon

%% Inital conditions for Simulation

global pathindex
pathindex = 1;
X0 = Track.X(1);
Y0 = Track.Y(1);
Psi0= -pi/2;
v0 = 1; %m/s
tp = 0.55;

vref = 1.25;


%% Run simulation and plot trajectory
pathindex = 1;

out = sim('CAR_SIM_MPC.slx',100);

%% Plots

close all
figure
plot(Track.X,Track.Y)
hold on
plot(out.X,out.Y)
axis equal
grid on
legend('Desired Trajectory', 'Vehicle Trajectory')
xlabel('X [m]')
ylabel('Y [m]')

figure
subplot(2,1,1)
plot(out.time,out.YawRateRef,out.time,out.YawRate)
grid on
legend('Reference Yaw Rate', 'Vehicle Yaw Rate')
xlabel('Time [s]')
ylabel('Yaw Rate [rad/s]')

subplot(2,1,2)
plot(out.time,out.deltaK)
grid on
axis tight
legend('Steering Command')
xlabel('Time [s]')
ylabel('Steering Angle [deg]')

figure
plot(out.time,out.SolverTime)
grid on
axis tight
xlabel('CHRONOS Compute Time [s]')
ylabel('Time [s]')


