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

%% Create MPC objects

% MIMO mpc
N = 10;         % Prediction Horizon
N_h_ctr = 5;    % Control Horizon
mpc = init_mpc(N,N_h_ctr);

% Horizontal Fan mpc
mpc_h = init_mpc(N_h_ctr);

% Vertical Fan mpc
mpc_v = init_mpc(N_h_ctr);

%% LPV system

% Get LPV model frozen at current state vector
[A,B,Bd,Ah,Bh,Av,Bv] = qLPV_TRMS_cascade_mpc_SS(Wh,Omh,Thth,Wv,Thtv);

% Initialize MIMO mpc system dynamics
mpc = init_mpc_system(mpc,eye(4)+Ts*A,Ts*B,Ts*Bd,eye(4),0,0);

% Initialize Horizontal Fan mpc system dynamics
mpc_h = init_mpc_system(mpc_h,1+Ts*Ah,Ts*Bh,0,1,0,0);

% Initialize Vorizontal Fan mpc system dynamics
mpc_v = init_mpc_system(mpc_v,1+Ts*Av,Ts*Bv,0,1,0,0);

%% Constraints

% Basic mpc constraints for reference
% x_min = [-2.9;-1;-1.7;-1.6;-0.6;-0.5];
% x_max = [2.9;1;1.2;1.6;0.6;1];
% u_min = [-2.5;-2];
% u_max = [2.5;2];

% MIMO mpc state constraints
x_min = [-1; -1.7; -0.6; -0.5];
x_max = [ 1;  1.2;  0.6;  1];
mpc = init_mpc_state_cnstr(mpc,x_min,x_max);

% MIMO mpc input constraints
u_min = [-2.9; -1.6];
u_max = [ 2.9;  1.6];
mpc = init_mpc_u_cnstr(mpc,u_min,u_max);

% Horizontal Fan state constraints
x_min = -2.9;
x_max = 2.9;
mpc_h = init_mpc_state_cnstr(mpc_h,x_min,x_max);

% Horizontal Fan input constraints
u_min = -2.5;
u_max = 2.5;
mpc_h = init_mpc_u_cnstr(mpc_h,u_min,u_max);

% Vertical Fan state constraints
x_min = -1.6;
x_max = 1.6;
mpc_v = init_mpc_state_cnstr(mpc_v,x_min,x_max);

% Vertical Fan input constraints
u_min = -2;
u_max = 2;
mpc_v = init_mpc_u_cnstr(mpc_v,u_min,u_max);

%% Terminal Ingredients

% MIMO mpc terminal ingredients computed using the dLQR method
Qx = diag([1 50 1 50]);     % State Penalty
Ru = 0.1;                   % Control Penalty
ter_constraint = 0;         % Only terminal cost
x_ref_is_y = 1;             % The terminal reference can be extracted 
                            % mpc tracking reference
mpc = init_mpc_ter_ingredients_dlqr(mpc,Qx,Ru,ter_constraint,x_ref_is_y);

% Horizontal Fan mpc terminal ingredients computed using the dLQR method
Qx = 50;                    % State Penalty
Ru = 1;                     % Control Penalty
ter_constraint = 0;         % Only terminal cost
x_ref_is_y = 0;             % The terminal reference cannnot be extracted 
                            % from the mpc tracking reference
mpc_h = init_mpc_ter_ingredients_dlqr(mpc_h,Qx,Ru,ter_constraint,x_ref_is_y);

% Vertical Fan mpc terminal ingredients computed using the dLQR method
Qx = 50;                    % State Penalty
Ru = 1;                     % Control Penalty
ter_constraint = 0;         % Only terminal cost
x_ref_is_y = 0;             % The terminal reference cannnot be extracted 
                            % from the mpc tracking reference
mpc_v = init_mpc_ter_ingredients_dlqr(mpc_v,Qx,Ru,ter_constraint,x_ref_is_y);

%% Costs

% MIMO mpc tracking penalty
Qe = diag([50 1 50 1]);
mpc = init_mpc_Tracking_cost(mpc,Qe);

% MIMO mpc control inputs variation penalty
Rdu = diag([5 5]);
mpc = init_mpc_DiffControl_cost(mpc,Rdu);

% Horizontal Fan mpc tracking penalty
Qe = 50;
mpc_h = init_mpc_Tracking_cost(mpc_h,Qe);

% Horizontal Fan mpc control inputs variation penalty
Rdu = 1;
mpc_h = init_mpc_DiffControl_cost(mpc_h,Rdu);

% Vertical Fan mpc tracking penalty
Qe = 50;
mpc_v = init_mpc_Tracking_cost(mpc_v,Qe);

% Vertical Fan mpc control inputs variation penalty
Rdu = 1;
mpc_v = init_mpc_DiffControl_cost(mpc_v,Rdu);


%% Init conditions for simulation

% Initialize optimization vector as all 0  
x0 = zeros(mpc.Nu+mpc.Nx,1);
u_prev = [0;0];     % init value for control
mpc.t = 500;        % increase t value to give preference to objectives
                    % over constraints

% Initialize optimization vector as all 0  
x0_h = zeros(mpc_h.Nu+mpc_h.Nx,1);
uh = 0;             % init value for control

% Initialize optimization vector as all 0  
x0_v = zeros(mpc_v.Nu+mpc_v.Nx,1);
uv = 0;             % init value for control
