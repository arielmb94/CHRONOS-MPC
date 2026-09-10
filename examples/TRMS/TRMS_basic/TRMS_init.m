%% Parameters

% Equilibrium angle for vertical dynamics
Thtv0 = -0.619426178368110;

% Initial state vector
x = [0.01,0.010,0.01,0.01,0.01,Thtv0+0.01]';

% Assign state vector variables    
Wh   = x(1);    % Horizontal Fan Angular Speed
Omh  = x(2);    % Horizontal Angular Rate
Thth = x(3);    % Horizontal Angle
Wv   = x(4);    % Vertical Fan Angular Speed
Omv  = x(5);    % Vertical Angular Rate
Thtv = x(6);    % Vertical Angle

% Sampling time
Ts = 0.1;
%% Create MPC object

N = 15;          % Prediction horizon: N*Ts = 1.5 s

% Create mpc struct
mpc = init_mpc(N);
%% Nominal LPV prediction model

% Freeze the LPV model at the initial operating point
sys = qLPV_TRMS_SS(Wh,Omh,Thth,Wv,Thtv);

% Initialize system dynamics
% System discretized with forward Euler discretization:
% x+ = (I+Ts*A)*x+Ts*B*u+Ts*Bd*d
mpc = init_mpc_dynamics(mpc,eye(6)+Ts*sys.A,Ts*sys.B,[]);

% All six states are tracked, so no separate output model is needed

%% Constraints

% State constraints
x_min = [-2.9;-1;-1.7;-1.6;-0.6;-0.5];
x_max = [2.9;1;1.2;1.6;0.6;1];
mpc = init_mpc_state_cnstr(mpc,x_min,x_max);

% Control input constraints
u_min = [-2.5;-2];
u_max = -u_min;
mpc = init_mpc_control_cnstr(mpc,u_min,u_max);

% Control inputs variation constraints
du_min = 0.2*u_min;
du_max = 0.2*u_max;
mpc = init_mpc_control_rate_cnstr(mpc,du_min,du_max);

%% Terminal ingredients

% DLQR weights used to compute the terminal cost matrix
Qx_dlqr = diag([1 1 50 1 1 1000]);
Ru_dlqr = 1;
x_ref_is_y = 1;

% The terminal matrix is computed from this nominal model
mpc = init_mpc_ter_ingredients_dlqr(mpc,Qx_dlqr,Ru_dlqr,x_ref_is_y);
%% Costs

% Tracking penalty
Qe = diag([1 50 1 1 500 1]);
mpc = init_mpc_Tracking_cost(mpc,Qe);

% Control inputs variation penalty
Rdu = diag([5 5]);
mpc = init_mpc_ControlRate_cost(mpc,Rdu);

% Control penalty
Ru = diag([1/2.5^2 1/2^2]);
mpc = init_mpc_Control_cost(mpc,Ru);

%% Finalize controller

u_prev = [0;0];
x_mpc = [Wh;Omh;Thth;Wv;Omv;Thtv-Thtv0];
mpc = build_chronos_mpc(mpc,x_mpc,u_prev);
