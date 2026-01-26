%% Car Parameters

m=1.1934; % mass of the vehicle
l_r=0.1049; % distance from center of gravity to rear wheel
L=0.174;
l_f=L-l_r; 

%% Create MPC object

N = 15;         % Prediction Horizon
N_h_ctr = 3;    % Control Horizon

% Create mpc struct
mpc = init_mpc(N,N_h_ctr);
%% LTI system

% Init speed
vx = 1;

% Compute Model to current vehicle speed
[A,B] = update_BM(vx);
C = [0 1];

% Initialize system dynamics
% System discretized with forward Euler discretization:
% x+ = (I+Ts*A)*x+Ts*B*u+Ts*Bd*d
mpc = init_mpc_system(mpc,eye(2)+Ts*A,Ts*B,0,C,0,0);

%% Constraints

% Control input constraints: +-40 degrees steering wheel angle
u_min = -0.7*ones(mpc.nu,1);
u_max = 0.7*ones(mpc.nu,1);
mpc = init_mpc_u_cnstr(mpc,u_min,u_max);

%% Costs

% Tracking penalty with Normalization
max_error = 3; % yaw rate max error [rad/s]

Qe = diag(100*ones(mpc.ny,1)/(max_error^2));
mpc = init_mpc_Tracking_cost(mpc,Qe);

% Control penalty with Normalization
max_u = 0.7; % steering wheel max anlge [rad]

Ru = 10/(max_u^2);    % Quadratic penalty on control action u'*Ru*u
ru = [];    % Linear penalty on control action vector: ru'*u 
mpc = init_mpc_Control_cost(mpc,Ru,ru);

%% Init conditions for simulation

% use warm start function to get optimization vector initial value
x_prev = [0; 0];
u_prev = 0;
x0 = init_mpc_warm_start(mpc,x_prev,u_prev);