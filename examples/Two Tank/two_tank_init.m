%% Parameters

% Tank area and gravity value
Ab = 1;
g = 9.81;

% Initial condition
h1 = 0.45;
h2 = 0.45;

% Sampling time
Ts = 0.01;

%% Create MPC object

N = 10;         % Prediction horizon: N*Ts = 0.1 s

% Create mpc struct
mpc = init_mpc(N);

% Limit optimization iterations per control sample
mpc.max_iter = 3;

%% Nominal prediction model

% Freeze the LPV model at the initial operating point
A = [-sqrt(2*g)*sqrt(h1)/(Ab*h1) 0;
     sqrt(2*g)*sqrt(h1)/(Ab*h1) -sqrt(2*g)*sqrt(h2)/(Ab*h2)];
B = [1/Ab; 0];

% Initialize system dynamics
% System discretized with forward Euler discretization:
% x+ = (I+Ts*A)*x+Ts*B*u+Ts*Bd*d
mpc = init_mpc_dynamics(mpc,eye(2)+Ts*A,Ts*B,[]);

% Replace the default full-state output: track only the second tank height
C = [0 1];

mpc = init_mpc_output(mpc,C,[],[]);

%% Constraints

% State constraints
x_min = 0.01*ones(mpc.nx,1);
x_max = 1*ones(mpc.nx,1);
mpc = init_mpc_state_cnstr(mpc,x_min,x_max);

% Control input constraints
u_min = 0*ones(mpc.nu,1);
u_max = 10*ones(mpc.nu,1);
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

%% Terminal Ingredients

% DLQR design weights used to compute the terminal cost matrix P
Qx_dlqr = diag([30 30]);
Ru_dlqr = 1;

% The tracked output is scalar, so pass the full-state terminal reference
x_ref_is_y = 0;

% P is computed once from the nominal model and is not updated online.
% When using a terminal cost, choose a representative initialization model.
mpc = init_mpc_ter_ingredients_dlqr(mpc,Qx_dlqr,Ru_dlqr,x_ref_is_y);

%% Costs

% Tracking penalty
Qe = diag(50*ones(mpc.ny,1));
mpc = init_mpc_Tracking_cost(mpc,Qe);

% Control inputs variation penalty
Rdu = 1;
mpc = init_mpc_ControlRate_cost(mpc,Rdu);

% Optional control-action cost (inactive in this example)
Ru = [];    % Quadratic penalty on control action u'*Ru*u
ru = [];    % Linear penalty on control action vector: ru'*u 
%mpc = init_mpc_Control_cost(mpc,Ru,ru);

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

%% Init conditions for simulation

% Finalize the controller for the initial state and previous input
x_prev = [h1; h2];
u_prev = 3.7;
mpc = build_chronos_mpc(mpc,x_prev,u_prev);
