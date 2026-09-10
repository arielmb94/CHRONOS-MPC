%% Parameters

% Tank reactor parameters
theta_f = 20;
k = 300;
M = 5;
xf = 0.3947;
xc = 0.3816;
alpha = 0.117;

% Initial state values
c0 = 0.2632;
v0 = 0.6519;
% State vector
x_prev = [c0;v0];

% Sampling Time
Ts = 0.1;

%% Create MPC object

N = 15;         % Prediction horizon: N*Ts = 1.5 s

% Create mpc struct
mpc = init_mpc(N);

%% Nominal prediction model

% Freeze the hybrid LPV model at the initial operating point
A = [-1/theta_f-k*exp(-M/v0) -k*c0*M*exp(-M/v0)/(v0^2);
     k*exp(-M/v0) -1/theta_f];
B = [0; -alpha*(v0-xc)];
Bd = [1/theta_f k*c0*M*exp(-M/v0)/(v0^2); xf/theta_f 0];

% Initialize system dynamics
% System discretized with forward Euler discretization:
% x+ = (I+Ts*A)*x+Ts*B*u+Ts*Bd*d
mpc = init_mpc_dynamics(mpc,eye(2)+Ts*A,Ts*B,Ts*Bd);

% CHRONOS tracks the full state by default, so no output model is needed

%% Constraints

% State constraints
x_min = [0;0];
x_max = [1;1];
slack_cost = 10;
mpc = init_mpc_state_cnstr(mpc,x_min,x_max,slack_cost,slack_cost);

% Control input constraints
u_min = 0;
u_max = 1;
mpc = init_mpc_control_cnstr(mpc,u_min,u_max);

% Control inputs variation constraints
du_min = -0.1*ones(mpc.nu,1);
du_max = 0.1*ones(mpc.nu,1);
mpc = init_mpc_control_rate_cnstr(mpc,du_min,du_max);

% Optional output constraints (inactive in this example)
y_min = [];
y_max = [];
%mpc = init_mpc_output_cnstr(mpc,y_min,y_max);

%% Optional custom constraints (inactive)

% Custom constrained signal: h = Ch*x+Dh*u+Dsuh*su+Ddh*dh
Ch = [];
Dh = [];
Dsuh = [];
Ddh = [];

h_min = [];
h_max = [];

%mpc = init_mpc_custom_cnstr(mpc,Ch,Dh,Dsuh,Ddh,h_min,h_max);

%% Optional terminal ingredients (inactive)

Qx_dlqr = [];
Ru_dlqr = [];

% With the default y = s, the tracking reference can also be used at the
% terminal stage when terminal ingredients are enabled
x_ref_is_y = 1;

%mpc = init_mpc_ter_ingredients_dlqr(mpc,Qx_dlqr,Ru_dlqr,x_ref_is_y);

%% Costs

% Tracking penalty
Qe = diag([5000 250]);
mpc = init_mpc_Tracking_cost(mpc,Qe);

% Optional control-rate cost (inactive; rate constraints remain active)
Rdu = [];
%mpc = init_mpc_ControlRate_cost(mpc,Rdu);

% Optional control-action cost (inactive)
Ru = [];
%mpc = init_mpc_Control_cost(mpc,Ru);

%% Optional custom cost (inactive)

% Custom performance signal: z = Cz*x+Dz*u+Dsuz*su+Ddz*dz
Cz = [];
Dz = [];
Dsuz = [];
Ddz = [];

Qz = [];    % Quadratic penalty on performance vector: z'*Qz*z
qz = [];    % Linear penalty on performance vector: qz'*z 

% Init performance cost
%mpc = init_mpc_Custom_cost(mpc,Cz,Dz,Dsuz,Ddz,Qz,qz);

%% Solver tuning

mpc.t = 500; % Barrier parameter tuned for this example (default: 50)

%% Finalize controller

u_prev = 0.45;
d = [1;v0];  % known input for the nominal model
mpc = build_chronos_mpc(mpc,x_prev,u_prev,d);
