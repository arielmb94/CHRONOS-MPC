%% Parameters

Ab = 1;
g = 9.81;

h1 = 0.45;
h2 = 0.45;

Ts = 0.01;

%% Create MPC object

N = 10;
N_h_ctr = 5;

mpc = init_mpc(N,N_h_ctr);
%% LTI system

A = [-sqrt(2*g)*sqrt(h1)/(Ab*h1) 0;
     sqrt(2*g)*sqrt(h1)/(Ab*h1) -sqrt(2*g)*sqrt(h2)/(Ab*h2)];

B = [1/Ab; 0];

C = [0 1];
%C = eye(2);

sys_ct = ss(A,B,C,0);
sys = c2d(sys_ct,Ts)

mpc = init_mpc_system(mpc,sys.A,sys.B,[],sys.C,sys.D,[]);

%% Constraints

x_min = 0.01*ones(mpc.nx,1);
x_max = 1*ones(mpc.nx,1);
mpc = init_mpc_state_cnstr(mpc,x_min,x_max);

x_ter_min = [];%;0.05*ones(nx,1);
x_ter_max = [];%1*ones(nx,1);
%mpc = init_mpc_ter_state_cnstr(mpc,x_ter_min,x_ter_max);

u_min = 0*ones(mpc.nu,1);
u_max = 10*ones(mpc.nu,1);
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

Qx = diag([30 30]);
Ru = 1;
ter_constraint = 0;
x_ref_is_y = 0;

[mpc] = init_mpc_ter_ingredients_dlqr(mpc,Qx,Ru,ter_constraint,x_ref_is_y);

%% Costs

Qe = diag(50*ones(mpc.ny,1));
mpc = init_mpc_Tracking_cost(mpc,Qe);

Rdu = 1;
mpc = init_mpc_DiffControl_cost(mpc,Rdu);

Ru = 1;
%mpc = init_mpc_Control_cost(mpc,Ru);

%% Performance Cost Matrix

Cz = [];
Dz = [];
Ddz = [];

Qz = [];

%mpc = init_mpc_LinPerf_cost(mpc,Qz,Cz,Dz,Ddz);

%%
%mpc = defLtiMpc(N,A,B,C,D,Bd,Dd,Qe,Rdu,Ru,Cz,Dz,Ddz,Qz,Ci,Di,Ddi,...
%    x_min,x_max,x_ter_min,x_ter_max,u_min,u_max,du_min,du_max,y_min,y_max,yi_min,yi_max)

%% Init conditions for simulation

x_prev = [h1; h2];
u_prev = 3.7;
x0 = init_mpc_warm_start(mpc,x_prev,u_prev);