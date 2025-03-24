%% Parameters

theta_f = 20;
k = 300;
M = 5;
xf = 0.3947;
xc = 0.3816;
alpha = 0.117;

Ts = 0.1;

c0 = 0.2632;
v0 = 0.6519;

x_prev = [c0;v0];

%% Create MPC object

N = 15;
N_h_ctr = 5;

mpc = init_mpc(N,N_h_ctr);
%% LTI system

A = [-1/theta_f-k*exp(-M/v0) -k*c0*M*exp(-M/v0)/(v0^2);
     k*exp(-M/v0) -1/theta_f];

B = [0; -alpha*(v0-xc)];

Bd = [1/theta_f k*c0*M*exp(-M/v0)/(v0^2); xf/theta_f 0];

C = [1 0];
C = eye(2);

sys_ct = ss(A,B,C,0);
sys = c2d(sys_ct,Ts)

mpc = init_mpc_system(mpc,eye(2)+Ts*A,Ts*B,Ts*Bd,C,[0;0],[]);

%% Constraints

x_min = [0;0];
x_max = [1;1];
mpc = init_mpc_state_cnstr(mpc,x_min,x_max);

x_ter_min = [];%;0.05*ones(nx,1);
x_ter_max = [];%1*ones(nx,1);
%mpc = init_mpc_ter_state_cnstr(mpc,x_ter_min,x_ter_max);

u_min = 0;
u_max = 1;
mpc = init_mpc_u_cnstr(mpc,u_min,u_max);

du_min = -0.1*ones(mpc.nu,1);
du_max = 0.1*ones(mpc.nu,1);
mpc = init_mpc_delta_u_cnstr(mpc,du_min,du_max);

y_min = [];
y_max = [];
%mpc = init_mpc_output_cnstr(mpc,y_min,y_max);

%% General Linear Inequalities

Ci = [];
Di = [];
Ddi = [];

yi_min = [];
yi_max = [];

%mpc = init_mpc_general_lin_ineq_cnstr(mpc,yi_min,yi_max,Ci,Di,Ddi);

%% Terminal Ingredients

Qx = diag([3000 0.1]);
Ru = 1;
ter_constraint = 0;
x_ref_is_y = 1;

%[mpc] = init_mpc_ter_ingredients_dlqr(mpc,Qx,Ru,ter_constraint,x_ref_is_y);

%% Costs

Qe = diag([5000 250]);
mpc = init_mpc_Tracking_cost(mpc,Qe);

Rdu = 1;
%mpc = init_mpc_DiffControl_cost(mpc,Rdu);

Ru = 10;
%mpc = init_mpc_Control_cost(mpc,Ru);

%% Performance Cost Matrix

Cz = [];
Dz = [];
Ddz = [];

Qz = [];

%mpc = init_mpc_LinPerf_cost(mpc,Qz,Cz,Dz,Ddz);
%% Set QP solver to be more aggressive

mpc.t = 500;

%%
%mpc = defLtiMpc(N,A,B,C,D,Bd,Dd,Qe,Rdu,Ru,Cz,Dz,Ddz,Qz,Ci,Di,Ddi,...
%    x_min,x_max,x_ter_min,x_ter_max,u_min,u_max,du_min,du_max,y_min,y_max,yi_min,yi_max)

%% Init conditions for simulation

x0 = 0.45*ones(mpc.Nu+mpc.Nx,1);
u_prev = 0.45;
d = [1;v0];
x0 = init_mpc_warm_start(mpc,x_prev,u_prev,d);
