%% Parameters
% Initial state vector: x = [theta1; theta2; dtheta1; dtheta2]
x = [0; pi; 0; 0];

% Sampling time
Ts = 0.02;

%% Create MPC object
N = 40;         % Prediction horizon: N*Ts = 0.8 s
mpc = init_mpc(N);
mpc.max_iter = 10;

%% Nominal LPV prediction model with Mass Matrix M
% Freeze the continuous LPV model and inertia matrix at the initial state
[Ac, Bc, Mc] = qLPV_Furuta_SS(x);

% Initialize system dynamics passing Mc directly to the discretization API
% choose between: _forward, _backward, or _tustin
mpc = init_mpc_dyn_continuous_forward(mpc, Ac, Bc, [], Ts, Mc);

C = eye(4);
mpc = init_mpc_output(mpc, C, [], []);

%% Constraints
% State constraints (Soft constraints with slack weight = 1e3)
x_min = [-3.1416; -100; -30; -30];
x_max = [ 3.1416;  100;  30;  30];
slack_cost = 1e3;
mpc = init_mpc_state_cnstr(mpc, x_min, x_max, slack_cost, slack_cost);

% Control input constraints: u in [-0.6, 0.6] N.m
u_min = -0.6;
u_max =  0.6;
mpc = init_mpc_control_cnstr(mpc, u_min, u_max);

%% Terminal ingredients
[Ac_eq, Bc_eq, Mc_eq] = qLPV_Furuta_SS([0; 0; 0; 0]);
Ad_eq = eye(4) + Ts * solveLinearSystemLU(Mc_eq, Ac_eq);
Bd_eq = Ts * solveLinearSystemLU(Mc_eq, Bc_eq);

Qx_dlqr = diag([10, 100, 5, 1]);
Ru_dlqr = 100;
x_ref_is_y = 1;

mpc = init_mpc_dynamics(mpc, Ad_eq, Bd_eq, []);
mpc = init_mpc_ter_ingredients_dlqr(mpc, Qx_dlqr, Ru_dlqr, x_ref_is_y);
mpc = init_mpc_dyn_continuous_forward(mpc, Ac, Bc, [], Ts, Mc);

%% Costs
Qe = diag([10, 100, 5, 1]);
mpc = init_mpc_Tracking_cost(mpc, Qe);

Ru = 100;
mpc = init_mpc_Control_cost(mpc, Ru);

%% build
u_prev = 0;
mpc = build_chronos_mpc(mpc, x, u_prev);