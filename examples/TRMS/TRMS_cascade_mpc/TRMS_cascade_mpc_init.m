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

%% Create MPC objects

N = 15;         % Prediction horizon: N*Ts = 1.5 s

% Outer MIMO MPC
mpc = init_mpc(N);

% Tail-rotor SISO MPC
mpc_h = init_mpc(N);

% Main-rotor SISO MPC
mpc_v = init_mpc(N);

% Limit the inner Newton iterations to control their online computation time
inner_max_iter = 4;
mpc_h.max_iter = inner_max_iter;
mpc_v.max_iter = inner_max_iter;

%% Nominal LPV prediction models

% Freeze the three LPV models at the initial operating point
[A,B,Bd,Ah,Bh,Av,Bv] = qLPV_TRMS_cascade_mpc_SS(Wh,Omh,Thth,Wv,Thtv);

% Outer body model: rotor references are the inputs and the main-rotor
% voltage is a measured disturbance
% System discretized with forward Euler discretization:
% x+ = (I+Ts*A)*x+Ts*B*u+Ts*Bd*d
mpc = init_mpc_dynamics(mpc,eye(4)+Ts*A,Ts*B,Ts*Bd);

% Tail-rotor model
mpc_h = init_mpc_dynamics(mpc_h,1+Ts*Ah,Ts*Bh,[]);

% Main-rotor model
mpc_v = init_mpc_dynamics(mpc_v,1+Ts*Av,Ts*Bv,[]);

%% Constraints

% Outer MIMO MPC: body-state and rotor-reference bounds
x_min_outer = [-1;-1.7;-0.6;-0.5];
x_max_outer = [1;1.2;0.6;1];
omega_ref_min = [-2.9;-1.6];
omega_ref_max = [2.9;1.6];
state_slack_cost = 10;
mpc = init_mpc_state_cnstr( ...
    mpc,x_min_outer,x_max_outer,state_slack_cost,state_slack_cost);
mpc = init_mpc_control_cnstr(mpc,omega_ref_min,omega_ref_max);

% Tail-rotor MPC: rotor-speed and motor-voltage bounds
Wh_min = -2.9;
Wh_max = 2.9;
uh_min = -2.5;
uh_max = 2.5;
mpc_h = init_mpc_state_cnstr( ...
    mpc_h,Wh_min,Wh_max,state_slack_cost,state_slack_cost);
mpc_h = init_mpc_control_cnstr(mpc_h,uh_min,uh_max);

% Main-rotor MPC: rotor-speed and motor-voltage bounds
Wv_min = -1.6;
Wv_max = 1.6;
uv_min = -2;
uv_max = 2;
mpc_v = init_mpc_state_cnstr( ...
    mpc_v,Wv_min,Wv_max,state_slack_cost,state_slack_cost);
mpc_v = init_mpc_control_cnstr(mpc_v,uv_min,uv_max);

%% Terminal ingredients

% The terminal matrices are computed once from these nominal models
Qx_outer_dlqr = diag([1 50 1 1000]);
Ru_outer_dlqr = 0.1;
outer_x_ref_is_y = 1;
mpc = init_mpc_ter_ingredients_dlqr( ...
    mpc,Qx_outer_dlqr,Ru_outer_dlqr,outer_x_ref_is_y);

% The inner MPCs receive their terminal references separately
Qx_h_dlqr = 50;
Ru_h_dlqr = 1;
Qx_v_dlqr = 50;
Ru_v_dlqr = 1;
inner_x_ref_is_y = 0;
mpc_h = init_mpc_ter_ingredients_dlqr( ...
    mpc_h,Qx_h_dlqr,Ru_h_dlqr,inner_x_ref_is_y);
mpc_v = init_mpc_ter_ingredients_dlqr( ...
    mpc_v,Qx_v_dlqr,Ru_v_dlqr,inner_x_ref_is_y);

%% Costs

% Outer MIMO MPC costs
Qe_outer = diag([50 1 500 1]);
Rdu_outer = diag([5 5]);
mpc = init_mpc_Tracking_cost(mpc,Qe_outer);
mpc = init_mpc_ControlRate_cost(mpc,Rdu_outer);

% Tail-rotor MPC costs
Qe_h = 50;
Rdu_h = 50;
mpc_h = init_mpc_Tracking_cost(mpc_h,Qe_h);
mpc_h = init_mpc_ControlRate_cost(mpc_h,Rdu_h);

% Main-rotor MPC costs
Qe_v = 50;
Rdu_v = 20;
mpc_v = init_mpc_Tracking_cost(mpc_v,Qe_v);
mpc_v = init_mpc_ControlRate_cost(mpc_v,Rdu_v);

%% Finalize controllers

omega_ref_prev = [0;0];
uh_prev = 0;
uv_prev = 0;

x_outer = [Omh;Thth;Omv;Thtv-Thtv0];
mpc = build_chronos_mpc(mpc,x_outer,omega_ref_prev,uv_prev);
outer_barrier_parameter = 500;
mpc.t = outer_barrier_parameter;

mpc_h = build_chronos_mpc(mpc_h,Wh,uh_prev);
mpc_v = build_chronos_mpc(mpc_v,Wv,uv_prev);
