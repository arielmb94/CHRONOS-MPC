%% Parameters

% Equlibrium angle for vertical dynamics
Thtv0 = -0.619426178368110;

% Initial state vector
x = [0.01,0.010,0.01,0.01,0.01,Thtv0+0.01]';

% Assign state vector variables    
Wh   = x(1);    % Horizontal Fan Angular Speed
Omh  = x(2);    % Horizontal Angular Rate
Thth = x(3);    % Horizontal Angle
Wv   = x(4);    % Vertical Fan Angular Speed
Thtv = x(6);    % Vertical Angle

% Sampling time
Ts = 0.1;
%% Create MPC object

N = 5;          % Prediction Horizon
N_h_ctr = 3;    % Control Horizon

% Create mpc struct
mpc = init_mpc(N,N_h_ctr);
%% LTI system

% Get LPV model frozen at current state vector
sys = qLPV_TRMS_SS(Wh,Omh,Thth,Wv,Thtv);
% Tracking for all 6 states
C = eye(6);

% Initialize system dynamics
% System discretized with forward Euler discretization:
% x+ = (I+Ts*A)*x+Ts*B*u+Ts*Bd*d
mpc = init_mpc_system(mpc,eye(6)+Ts*sys.A,Ts*sys.B,0,C,0,0);

%% Constraints

% State constraints
x_min = [-2.9;-1;-1.7;-1.6;-0.6;-0.5];
x_max = [2.9;1;1.2;1.6;0.6;1];
mpc = init_mpc_state_cnstr(mpc,x_min,x_max);

% State constraints only on terminal states
x_ter_min = [-2.9;-1;-1.7;-1.6;-0.6;-0.5];
x_ter_max = [2.9;1;1.2;1.6;0.6;1];
%mpc = init_mpc_ter_state_cnstr(mpc,x_ter_min,x_ter_max);

% Control input constraints
u_min = [-2.5;-2];
u_max = -u_min;
mpc = init_mpc_u_cnstr(mpc,u_min,u_max);

% Control inputs variation constraints
du_min = -0.5*ones(mpc.nu,1);
du_max = 0.5*ones(mpc.nu,1);
%mpc = init_mpc_delta_u_cnstr(mpc,du_min,du_max);

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
Qx = diag([1 50 1 1 50 1]);     % State Penalty
Ru = 1;                         % Control Penalty
ter_constraint = 0;             % Only terminal cost
x_ref_is_y = 1;                 % The terminal reference can be extracted 
                                % mpc tracking reference

% Initialize terminal ingredients using the dLQR method
mpc = init_mpc_ter_ingredients_dlqr(mpc,Qx,Ru,ter_constraint,x_ref_is_y);
%% Costs

% Tracking penalty
Qe = diag([1 50 1 1 50 1]);
mpc = init_mpc_Tracking_cost(mpc,Qe);

% Control inputs variation penalty
Rdu = diag([1 1]);
mpc = init_mpc_DiffControl_cost(mpc,Rdu);

% Control penalty
%Ru = diag([0.5 1]);
%mpc = init_mpc_Control_cost(mpc,Ru);

%% Init conditions for simulation

u_prev = [0;0];
% use warm start function to get optimization vector initial value
x0 = init_mpc_warm_start(mpc,x,u_prev);