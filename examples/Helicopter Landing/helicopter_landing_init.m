clear all
close all
clc

%% Parameters

kx = 0.05;
ky = 0.02;

g = -9.8;

Ts = 0.1;

initial_px = 44.7;   % Initial position in the X axis
initial_vx = 13.2;  % Initial velocity in X axis
initial_py = 40; % Initial position in the Y axis
initial_vy = -14.75; % Initial velocity in Y axis

x_prev = [initial_px; initial_vx; initial_py; initial_vy];

%% Create MPC object

N = 15;
N_h_ctr = 15;

mpc = init_mpc(N,N_h_ctr);
%% LTI system

rho_1 = -kx * initial_vx;
rho_2 = -ky * initial_vy;

A = [1 1*Ts 0 0; 
    0 1+Ts*rho_1 0 0;
    0 0 1 1*Ts;
    0 0 0 1+Ts*rho_2];

B = [0 0;
     Ts 0;
     0 0;
     0 Ts];

Bd = [0;
      0;
      0;
      Ts];

C = eye(4);

mpc = init_mpc_system(mpc, A, B, Bd, C, [0;0;0;0], []);

%% Constraints

x_min = [0;-20;0;-20];
x_max = [1000;20;1000;20];
mpc = init_mpc_state_cnstr(mpc,x_min,x_max);

x_ter_min = [];%;0.05*ones(nx,1);
x_ter_max = [];%1*ones(nx,1);
%mpc = init_mpc_ter_state_cnstr(mpc,x_ter_min,x_ter_max);

u_min = [-20; -20];
u_max = [20; 20];
mpc = init_mpc_u_cnstr(mpc,u_min,u_max);

% du_min = -0.1*ones(mpc.nu,1);
% du_max = 0.1*ones(mpc.nu,1);
% mpc = init_mpc_delta_u_cnstr(mpc,du_min,du_max);

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

Qx = diag([1, 0.3, 0.7, 0.3]);
Ru = [0.2, 0.2];
ter_constraint = 0;
x_ref_is_y = 1;

%[mpc] = init_mpc_ter_ingredients_dlqr(mpc,Qx,Ru,ter_constraint,x_ref_is_y);

%% Costs

Qe = diag([1, 0.3, 0.7, 0.3]);
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

%x0 = 0.45*ones(mpc.Nu+mpc.Nx,1);
u_prev = [0;0];
d = [g];
x0 = init_mpc_warm_start(mpc,x_prev,u_prev,d);
