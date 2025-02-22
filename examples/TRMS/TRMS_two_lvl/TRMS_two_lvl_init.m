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

% MIMO mpc
N = 10;
N_h_ctr = 5;

mpc = init_mpc(N,N_h_ctr);

% Horizontal Fan mpc
mpc_h = init_mpc(N_h_ctr);

% Vertical Fan mpc
mpc_v = init_mpc(N_h_ctr);
%% LTI system

[A,B,Bd,Ah,Bh,Av,Bv] = qLPV_TRMS_two_lvl_SS(Wh,Omh,Thth,Wv,Thtv);

mpc = init_mpc_system(mpc,eye(4)+Ts*A,Ts*B,Ts*Bd,eye(4),0,0);

mpc_h = init_mpc_system(mpc_h,1+Ts*Ah,Ts*Bh,0,1,0,0);

mpc_v = init_mpc_system(mpc_v,1+Ts*Av,Ts*Bv,0,1,0,0);

%% Constraints

% x_min = [-2.9;-1;-1.7;-1.6;-0.6;-0.5];
% x_max = [2.9;1;1.2;1.6;0.6;1];
% u_min = [-2.5;-2];
% u_max = [2.5;2];

x_min = [-1; -1.7; -0.6; -0.5];
x_max = [ 1;  1.2;  0.6;  1];
mpc = init_mpc_state_cnstr(mpc,x_min,x_max);

x_ter_min = [];%;0.05*ones(nx,1);
x_ter_max = [];%1*ones(nx,1);
%mpc = init_mpc_ter_state_cnstr(mpc,x_ter_min,x_ter_max);

u_min = [-2.9; -1.6];
u_max = [ 2.9;  1.6];
mpc = init_mpc_u_cnstr(mpc,u_min,u_max);

du_min = -0.5*ones(mpc.nu,1);
du_max = 0.5*ones(mpc.nu,1);
%mpc = init_mpc_delta_u_cnstr(mpc,du_min,du_max);

y_min = [];
y_max = [];
%mpc = init_mpc_output_cnstr(mpc,y_min,y_max);

x_min = -2.9;
x_max = 2.9;
mpc_h = init_mpc_state_cnstr(mpc_h,x_min,x_max);

u_min = -2.5;
u_max = 2.5;
mpc_h = init_mpc_u_cnstr(mpc_h,u_min,u_max);

x_min = -1.6;
x_max = 1.6;
mpc_v = init_mpc_state_cnstr(mpc_v,x_min,x_max);

u_min = -2;
u_max = 2;
mpc_v = init_mpc_u_cnstr(mpc_v,u_min,u_max);

%% General Linear Inequalities

Ci = [];
Di = [];
Ddi = [];

yi_min = [];
yi_max = [];

%mpc = init_mpc_general_lin_ineq_cnstr(mpc,yi_min,yi_max,Ci,Di,Ddi);

%% Terminal Ingredients

Qx = diag([1 50 1 50]);
Ru = 0.1;
ter_constraint = 0;
x_ref_is_y = 1;

mpc = init_mpc_ter_ingredients_dlqr(mpc,Qx,Ru,ter_constraint,x_ref_is_y);

Qx = 50;
Ru = 1;
ter_constraint = 0;
x_ref_is_y = 0;

mpc_h = init_mpc_ter_ingredients_dlqr(mpc_h,Qx,Ru,ter_constraint,x_ref_is_y);

Qx = 50;
Ru = 1;
ter_constraint = 0;
x_ref_is_y = 0;

mpc_v = init_mpc_ter_ingredients_dlqr(mpc_v,Qx,Ru,ter_constraint,x_ref_is_y);

%% Costs

Qe = diag([50 1 50 1]);
mpc = init_mpc_Tracking_cost(mpc,Qe);

Rdu = diag([5 5]);
mpc = init_mpc_DiffControl_cost(mpc,Rdu);


Qe = 50;
mpc_h = init_mpc_Tracking_cost(mpc_h,Qe);

Rdu = 1;
mpc_h = init_mpc_DiffControl_cost(mpc_h,Rdu);


Qe = 50;
mpc_v = init_mpc_Tracking_cost(mpc_v,Qe);

Rdu = 1;
mpc_v = init_mpc_DiffControl_cost(mpc_v,Rdu);

%Ru = diag([0.5 1]);
%mpc = init_mpc_Control_cost(mpc,Ru);

%% Performance Cost Matrix

Cz = [-1 0 0 0 0 0;  
      0 0 0 -1 0 0];
Dz = [0 0 1 0;
      0 0 0 1];
Ddz = 0;

Qz = diag([100 100]);

%mpc = init_mpc_LinPerf_cost(mpc,Cz,Dz,Ddz,Qz);

%% Init conditions for simulation

x0 = zeros(mpc.Nu+mpc.Nx,1);
u_prev = [0;0];

mpc.t = 500;

x0_h = zeros(mpc_h.Nu+mpc_h.Nx,1);
uh = 0;

mpc_h.t = 50;

x0_v = zeros(mpc_v.Nu+mpc_v.Nx,1);
uv = 0;

mpc_v.t = 50;
