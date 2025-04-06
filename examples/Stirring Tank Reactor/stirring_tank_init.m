%% Parameters

% Tank reactor parameters
theta_f = 20;
k = 300;
M = 5;
xf = 0.3947;
xc = 0.3816;
alpha = 0.117;

% Intial state values
c0 = 0.2632;
v0 = 0.6519;
% State vector
x_prev = [c0;v0];

% Sampling Time
Ts = 0.1;

%% Create MPC object

N = 15;         % Prediction Horizon
N_h_ctr = 5;    % Control Horizon

% Create mpc struct
mpc = init_mpc(N,N_h_ctr);

%% Define system dynamics

% Get LPV model frozen at current state vector
A = [-1/theta_f-k*exp(-M/v0) -k*c0*M*exp(-M/v0)/(v0^2);
     k*exp(-M/v0) -1/theta_f];
B = [0; -alpha*(v0-xc)];
Bd = [1/theta_f k*c0*M*exp(-M/v0)/(v0^2); xf/theta_f 0];
C = eye(2);

% Initialize system dynamics
mpc = init_mpc_system(mpc,eye(2)+Ts*A,Ts*B,Ts*Bd,C,[0;0],0);

%% Constraints

% State constraints
x_min = [0;0];
x_max = [1;1];
mpc = init_mpc_state_cnstr(mpc,x_min,x_max);

% State constraints only on terminal states
x_ter_min = [];
x_ter_max = [];
%mpc = init_mpc_ter_state_cnstr(mpc,x_ter_min,x_ter_max);

% Control input constraints
u_min = 0;
u_max = 1;
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
% yi = Ci*x+Di*u+Ddi*di
Ci = [];
Di = [];
Ddi = [];

yi_min = [];
yi_max = [];

%mpc = init_mpc_general_lin_ineq_cnstr(mpc,yi_min,yi_max,Ci,Di,Ddi);

%% Terminal Ingredients

% Terminal ingredients are computed using the dLQR method
Qx = [];                % State Penalty
Ru = [];                % Control Penalty
ter_constraint = 0;     % Only terminal cost
x_ref_is_y = 1;         % The terminal reference can be extracted 
                        % mpc tracking reference

% Initialize terminal ingredients using the dLQR method
%[mpc] = init_mpc_ter_ingredients_dlqr(mpc,Qx,Ru,ter_constraint,x_ref_is_y);

%% Costs

% Tracking penalty
Qe = diag([5000 250]);
mpc = init_mpc_Tracking_cost(mpc,Qe);

% Control inputs variation penalty
Rdu = [];
%mpc = init_mpc_DiffControl_cost(mpc,Rdu);

% Control penalty
Ru = [];
%mpc = init_mpc_Control_cost(mpc,Ru);

%% Performance Cost Matrix

% Performance Vector are defined as:
% z = Cz*x+Dz*u+Ddz*dz
Cz = [];
Dz = [];
Ddz = [];

Qz = [];    % Quadratic penalty on performance vector: z'*Qz*z
qz = [];    % Linear penalty on performance vector: qz'*z 

% Init performance cost
%mpc = init_mpc_LinPerf_cost(mpc,Cz,Dz,Ddz,Qz,qz);

%% Set QP solver to be more aggressive

mpc.t = 500; % Default value is t = 50, increasing t makes the solver give 
             % priority to cost penalties over constraints

%% Init conditions for simulation

% use warm start function to get optimization vector initial value
u_prev = 0.45;
d = [1;v0];
x0 = init_mpc_warm_start(mpc,x_prev,u_prev,d);
