%% Call the mpc problem initialization script
stirring_tank_init

%% Define simulation duration and reference parameters

% Duration
tsim = 240; % seconds
Sim_samples = tsim/Ts;
time = 0:Ts:tsim-Ts;

% Define concentration state reference
clear ref_c_vec
ref_c_vec = zeros(1,Sim_samples);
ref_c_vec(time < 90) = 0.27 + (0.65 - 0.27) * time(time < 90) / 90;
ref_c_vec(time >= 90 & time < 180) = 0.65;
ref_c_vec(time >= 180) = 0.65 - (0.65 - 0.3) * (time(time >= 180) - 180) / 60;

xref = [ref_c_vec(1);- M/(log(1/(theta_f*k*ref_c_vec(1))*(1-ref_c_vec(1))))];

%% Run Simulation

my_sched_fun = @(x) calc_rho_stirring_tank(x, k, M);
my_jacob_fun = @(x) calc_jacob_stirring_tank(x, k, M);

% clear storage variables
clear c_dat v_dat u_dat ref_v_dat t_dat

% Simulation Loop
for i = 1:Sim_samples

% Update concentration state reference
ref_c = ref_c_vec(i);
% Compute temperature state reference from steady state equilibrium 
% equation
ref_v = - M/(log(1/(theta_f*k*ref_c)*(1-ref_c)));

% Assign state vector variables  
ck = x_prev(1);
vk = x_prev(2);

% Number of schedulling variables
n_rho = 3;

tic
% Update LPV model
% System discretized with forward Euler discretization:
% x+ = (I+Ts*A)*x+Ts*B*u+Ts*Bd*d

% A = A0 + A1*rho_1 ...
A0 = eye(2) + Ts * [-1/theta_f, 0; 0, -1/theta_f];
A1 = Ts * [-1,  0; 1, 0];
A2 = Ts * [ 0, -1; 0, 0];
A3 = zeros(2,2);

% B = B0 + B1*rho_1 ...
B0 = Ts * [0; alpha*xc];
B1 = zeros(2,1);
B2 = zeros(2,1);
B3 = Ts * [0; -alpha];

% Bd = Bd0 + Bd1*rho_1 ...
Bd0 = Ts * [1/theta_f, 0; xf/theta_f, 0];
Bd1 = zeros(2,2);
Bd2 = Ts * [0, 1; 0, 0];
Bd3 = zeros(2,2);

n_iter = 5;
% [Pk, x0] = compute_schedul_iterative(mpc, x0, x_prev, u_prev, xref, d, ...
%                                               my_sched_fun, n_rho, n_iter, ...
%                                               A0, B0, Bd0, ...
%                                               A1, A2, A3, ...
%                                               B1, B2, B3, ...
%                                               Bd1, Bd2, Bd3);
%Pk = compute_schedul_frozen(mpc, x_prev, my_sched_fun);
Pk = compute_schedul_recursive(mpc, x0, x_prev, n_rho, ...
                               my_sched_fun, my_jacob_fun, ...
                               my_sched_fun(x_min), my_sched_fun(x_max));


% A_lpv = eye(2)+Ts*[-1/theta_f-k*exp(-M/vk) -k*ck*M*exp(-M/vk)/(vk^2);
%      k*exp(-M/vk) -1/theta_f];
% B_lpv = Ts*[0; -alpha*(vk-xc)];
% Bd_lpv = Ts*[1/theta_f k*ck*M*exp(-M/vk)/(vk^2); xf/theta_f 0];
% Update mpc problem dynamics
mpc = update_mpc_sys_dynamics(mpc, A0, B0, Bd0, Pk, n_rho, ...
                                  A1, A2, A3, ...      
                                  B1, B2, B3, ...       
                                  Bd1, Bd2, Bd3);

% Update mpc disturbance vector
d = [1;vk];
% Update reference vector
xref = [ref_c;ref_v];
% Solve mpc iteration
%[u_prev,x0] = mpc_solve(mpc,x0,x_prev,u_prev,xref,d,[],[],[]);
[u_prev,x0] = mpc_solve(mpc,x0,x_prev,u_prev,xref,d,[],[],[], ...
                        A0, Bd0, Pk, n_rho, ...
                        A1, A2, A3, ...
                        B1, B2, B3, ...
                        Bd1, Bd2, Bd3);
tk = toc;

% Store variables values for plotting and analysis
c_dat(:,i) = ck;
v_dat(:,i) = vk;
u_dat(i) = u_prev;
ref_v_dat(i) = ref_v;
t_dat(i) = tk;

% Forward Euler step of Stirring Tank nonlinear dynamics
ck = ck + Ts*((1-ck)/theta_f - k*ck*exp(-M/vk));
vk = vk + Ts*((xf-vk)/theta_f + k*ck*exp(-M/vk)-alpha*u_prev*(vk-xc));
% update state vector for the following iteration
x_prev = [ck;vk];

end

%% Plots
close all

figure
ax1 = subplot(3,1,1);
plot(time,ref_c_vec,'r',time,c_dat,'b')
grid on
ylim([0.2 0.7])
xlim([0 240])
legend('Concentration Ref.','c_k')
xlabel('Time (s)')
ylabel('Tank Concentration')

ax2 = subplot(3,1,2);
plot(time,ref_v_dat,'r',time,v_dat,'b')
grid on
ylim([0.5 0.7])
xlim([0 240])
legend('Temperature Ref.','v_k')
xlabel('Time (s)')
ylabel('Tank Temperature')

ax3 = subplot(3,1,3);
plot(time,u_dat,time(1:end-1),diff(u_dat))
grid on
xlim([0 240])
legend('u','\Delta u')
xlabel('Time (s)')
ylabel('Coolant Flow Rate')

linkaxes([ax1,ax2,ax3],'x')

figure
plot(time,t_dat)
xlim([0 240])
title('Compute Time (s)')
xlabel('Time (s)')
grid on