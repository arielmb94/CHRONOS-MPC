%% Load parameters and controllers

clearvars
close all

% Sampling time used throughout the simulator
Ts = 1/50; % s

% Identified rear-motor parameters
motor_data = load('xmotor.mat');
xmotor = motor_data.x;

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

% Identified SAV and tyre parameters
vehicle_data = load('xNL.mat');
vehicle_parameters = vehicle_data.x;

C_sigma_2 = vehicle_parameters(1);
C_sigma_1 = vehicle_parameters(2);
C_sigma_0 = vehicle_parameters(3);
C_alpha_f2 = vehicle_parameters(4);
C_alpha_f1 = vehicle_parameters(5);
C_alpha_f0 = vehicle_parameters(6);
C_alpha_r2 = vehicle_parameters(7);
C_alpha_r1 = vehicle_parameters(8);
C_alpha_r0 = vehicle_parameters(9);
eps = vehicle_parameters(10);
wn = vehicle_parameters(11);
tau = vehicle_parameters(12);
Iz = vehicle_parameters(13);
R = vehicle_parameters(14);

% Steering-actuator model, including the identified input delay
servoSIM = ss([0 1; -wn^2 -2*eps*wn],[0; wn^2],[1 0],0);
servoSIM.InputDelay = tau;

% Disturbance observer for the motor-control layer
Init_Observer_Motor_PI

% Existing longitudinal and motor controllers
load("controllers.mat")

% CHRONOS lateral MPC
sav_lat_controller_init

%% Load planned path

% Generate the Track structure used by the path-planning layer
Cartoon

%% Initial simulation conditions

global pathindex
pathindex = 1;
X0 = Track.X(1);
Y0 = Track.Y(1);
Psi0 = -pi/2; % rad
v0 = 1;       % Initial longitudinal speed (m/s)
tp = 0.55;    % Pure-Pursuit look-ahead time (s)
vref = 1.5;  % Longitudinal speed reference (m/s)

%% Run simulation

pathindex = 1;
simulation_time = 100; % s
out = sim('CAR_SIM_MPC.slx',simulation_time);

%% Plot trajectory and lateral-control signals

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
subplot(3,1,1)
plot(out.time,out.Vx)
grid on
title('Scheduling signal: vx')
xlabel('Time [s]')
ylabel('Vehcile Speed [m/s]')

subplot(3,1,2)
plot(out.time,out.YawRateRef,out.time,out.YawRate)
grid on
title('Yaw Rate Tracking')
legend('Reference Yaw Rate', 'Vehicle Yaw Rate')
xlabel('Time [s]')
ylabel('Yaw Rate [rad/s]')

subplot(3,1,3)
plot(out.time,out.deltaK)
grid on
axis tight
legend('Steering Command')
xlabel('Time [s]')
ylabel('Steering Angle [deg]')

figure
subplot(2,1,1)
plot(out.time,out.SolverTime)
grid on
axis tight
xlabel('Simulation Time [s]')
ylabel('Online Controller Time [s]')

subplot(2,1,2)
stairs(out.time,out.iter)
grid on
axis tight
xlabel('Simulation Time [s]')
ylabel('Newton Iterations')
