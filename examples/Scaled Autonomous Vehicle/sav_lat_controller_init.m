%% Vehicle parameters

m = 1.1934;     % Vehicle mass (kg)
l_r = 0.1049;   % Distance from center of mass to rear axle (m)
L = 0.174;      % Wheelbase (m)
l_f = L-l_r;    % Distance from center of mass to front axle (m)

%% Create MPC object

N = 50; % Prediction horizon: N*Ts = 1 s

% This template initializes the persistent MPC state used online
mpc_initial = init_mpc(N);
mpc_initial.max_iter = 3;
%% Initial LPV prediction model

% Freeze the bicycle model at the initial longitudinal speed
vx = 1; % m/s
[A,B] = update_BM(vx);

% Forward-Euler discretization:
% x+ = (I+Ts*A)*x+Ts*B*u+Ts*Bd*d
mpc_initial = init_mpc_dynamics(mpc_initial,eye(2)+Ts*A,Ts*B,[]);

% Track yaw rate only: y = C*x = yaw_rate
C_yaw_rate = [0 1];
mpc_initial = init_mpc_output(mpc_initial,C_yaw_rate,[],[]);

%% Constraints

% Steering-command bounds: approximately +/-40 degrees
steering_limit = 0.7; % rad
steering_min = -steering_limit;
steering_max = steering_limit;
mpc_initial = init_mpc_control_cnstr(mpc_initial,steering_min,steering_max);

%% Costs

% Yaw-rate tracking cost with normalization
yaw_rate_error_scale = 3; % rad/s
yaw_rate_tracking_weight = 100;

Qe = yaw_rate_tracking_weight/yaw_rate_error_scale^2;

mpc_initial = init_mpc_Tracking_cost(mpc_initial,Qe);

% Steering cost with normalization
steering_scale = steering_limit;
steering_weight = 10;
Ru = steering_weight/steering_scale^2;
linear_steering_weight = [];

mpc_initial = init_mpc_Control_cost(mpc_initial,Ru,linear_steering_weight);

%% Finalize controller

% Initial lateral state and previous steering command
x_prev = [0;0];
u_prev = 0;
mpc_initial = build_chronos_mpc(mpc_initial,x_prev,u_prev);
