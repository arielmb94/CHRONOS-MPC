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

N = 10;         % Prediction Horizon
N_h_ctr = 5;    % Control Horizon

% Create mpc struct
mpc = init_mpc(N,N_h_ctr);
%% LTI system

% Get LPV model frozen at current state vector
A = [-sqrt(2*g)*sqrt(h1)/(Ab*h1) 0;
     sqrt(2*g)*sqrt(h1)/(Ab*h1) -sqrt(2*g)*sqrt(h2)/(Ab*h2)];
B = [1/Ab; 0];

% Tracking objective is the water height on the second tank
C = [0 1];

% Initialize system dynamics
% System discretized with forward Euler discretization:
% x+ = (I+Ts*A)*x+Ts*B*u+Ts*Bd*d
mpc = init_mpc_system(mpc,eye(2)+Ts*A,Ts*B,0,C,0,0);

%% Constraints

% State constraints
x_min = 0.01*ones(mpc.nx,1);
x_max = 1*ones(mpc.nx,1);
mpc = init_mpc_state_cnstr(mpc,x_min,x_max);

% State constraints only on terminal states
x_ter_min = [];
x_ter_max = [];
%mpc = init_mpc_ter_state_cnstr(mpc,x_ter_min,x_ter_max);

% Control input constraints
u_min = 0*ones(mpc.nu,1);
u_max = 10*ones(mpc.nu,1);
mpc = init_mpc_u_cnstr(mpc,u_min,u_max);

% Control inputs variation constraints
du_min = -0.1*ones(mpc.nu,1);
du_max = 0.1*ones(mpc.nu,1);
mpc = init_mpc_delta_u_cnstr(mpc,du_min,du_max);

% Output constraints
y_min = [];
y_max = [];
%mpc = init_mpc_output_cnstr(mpc,y_min,y_max);

%% General Linear Inequalities

% General Linear Inequalities are defined as:
% yh = Ch*x+Dh*u+Ddh*di
Ch = [];
Dh = [];
Ddh = [];

h_min = [];
h_max = [];

%mpc = init_mpc_lin_custom_cnstr(mpc,h_min,h_max,Ch,Dh,Ddh);

%% Terminal Ingredients

% Terminal ingredients are computed using the dLQR method
Qx = diag([30 30]);         % State Penalty
Ru = 1;                     % Control Penalty
ter_constraint = 0;         % Only terminal cost
x_ref_is_y = 0;             % The terminal reference cannnot be extracted 
                            % from the mpc tracking reference

% Initialize terminal ingredients using the dLQR method                           
[mpc] = init_mpc_ter_ingredients_dlqr(mpc,Qx,Ru,ter_constraint,x_ref_is_y);

%% Costs

% Tracking penalty
Qe = diag(50*ones(mpc.ny,1));
mpc = init_mpc_Tracking_cost(mpc,Qe);

% Control inputs variation penalty
Rdu = 1;
mpc = init_mpc_DiffControl_cost(mpc,Rdu);

% Control penalty
Ru = [];    % Quadratic penalty on control action u'*Ru*u
ru = [];    % Linear penalty on control action vector: ru'*u 
%mpc = init_mpc_Control_cost(mpc,Ru,ru);

%% Performance Cost Matrix

% Performance Vector are defined as:
% z = Cz*x+Dz*u+Ddz*dz
Cz = [];
Dz = [];
Ddz = [];

Qz = [];    % Quadratic penalty on performance vector: z'*Qz*z
qz = [];    % Linear penalty on performance vector: qz'*z 

% Init performance cost
%mpc = init_mpc_Lin_Custom_cost(mpc,Cz,Dz,Ddz,Qz,qz);

%% Init conditions for simulation

% use warm start function to get optimization vector initial value
x_prev = [h1; h2];
u_prev = 3.7;
x0 = init_mpc_warm_start(mpc,x_prev,u_prev);