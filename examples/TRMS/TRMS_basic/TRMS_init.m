%% Parameters

Thtv0 = -0.619426178368110;
x = [0.01,0.010,0.01,0.01,0.01,Thtv0+0.01]';

Wh   = x(1);
Omh  = x(2);
Thth = x(3);
Wv   = x(4);
Thtv = x(6);

Ts = 0.1;

%% Create MPC object

N = 5;
N_h_ctr = 3;

mpc = init_mpc(N,N_h_ctr);
%% LTI system

sys = qLPV_TRMS_SS(Wh,Omh,Thth,Wv,Thtv);
C = eye(6);
% C = [1 0 0 0 0 0;
%     0 0 1 0 0 0;
%     0 0 0 1 0 0;
%     0 0 0 0 0 1];

mpc = init_mpc_system(mpc,eye(6)+Ts*sys.A,Ts*sys.B,0,C,0,0);

%% Constraints

x_min = [-2;-2;-2;-2;-0.6;-0.4];
x_max = -x_min;
mpc = init_mpc_state_cnstr(mpc,x_min,x_max);

x_ter_min = [];%;0.05*ones(nx,1);
x_ter_max = [];%1*ones(nx,1);
%mpc = init_mpc_ter_state_cnstr(mpc,x_ter_min,x_ter_max);

u_min = [-2.5;-2];
u_max = -u_min;
mpc = init_mpc_u_cnstr(mpc,u_min,u_max);

du_min = -0.5*ones(mpc.nu,1);
du_max = 0.5*ones(mpc.nu,1);
%mpc = init_mpc_delta_u_cnstr(mpc,du_min,du_max);

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

Qx = diag([1 50 1 1 50 1]);
Ru = 0.1;
ter_constraint = 0;
x_ref_is_y = 1;

[mpc] = init_mpc_ter_ingredients_dlqr(mpc,Qx,Ru,ter_constraint,x_ref_is_y);

%% Costs

Qe = diag([1 50 1 1 50 1]);
mpc = init_mpc_Tracking_cost(mpc,Qe);

Rdu = diag([1 1]);
mpc = init_mpc_DiffControl_cost(mpc,Rdu);

%Ru = diag([0.5 1]);
%mpc = init_mpc_Control_cost(mpc,Ru);

%% Init conditions for simulation

x0 = zeros(mpc.Nu+mpc.Nx,1);
u_prev = [0;0];

mpc.t = 500;
