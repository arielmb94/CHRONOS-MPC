%% Call the mpc problem initialization script

TRMS_init
%% Define simulation duration and reference parameters

% Duration
tsim = 200; % seconds
Sim_samples = tsim/Ts;
time = 0:Ts:tsim-Ts;

% Horizontal Angle sinousidal reference parameters
freq_TththRef = 1/31; % Hz
offset_TththRef = 0;
ampl_TththRef = 1;

% Vertical Angle sinousidal reference parameters
freq_TthtvRef = 1/47; % Hz
offset_TthtvRef = -0.6;
ampl_TthtvRef = 0.5;

% sinousidal references
TththRef_v = offset_TththRef + ampl_TththRef*sin(2*pi*(freq_TththRef)*time);
TthtvRef_v = offset_TthtvRef + ampl_TthtvRef*sin(2*pi*(freq_TthtvRef)*time);

%% LPV Setup for new methods
n = mpc.N;
n_rho = 5; % Wh, Omh, Thth, Wv, Thtv
n_iter = 3; % For SQP method

% Define how to extract rho from the MPC state vector 
% x_mpc = [Wh; Omh; Thth; Wv; Omv; Thtv - Thtv0]
% Note: We add Thtv0 back to state 6 to get the physical Thtv for the LPV function
my_sched_fun = @(x_st) [x_st(1); x_st(2); x_st(3); x_st(4); x_st(6) + Thtv0];

% Interface functions for traj_mat (Defined at the end of the script)
compute_A = @(rho) get_discrete_A(rho, Ts);
compute_B = @(rho) get_discrete_B(rho, Ts);

% Jacobian and bounds for Recursive method
my_jacob_fun = @(x) [1 0 0 0 0 0; 
                     0 1 0 0 0 0; 
                     0 0 1 0 0 0; 
                     0 0 0 1 0 0; 
                     0 0 0 0 0 1];
% Physical bounds based on x_min and x_max from TRMS_init
rho_min = [-2.9; -1.0; -1.7; -1.6; -0.5 + Thtv0];
rho_max = [ 2.9;  1.0;  1.2;  1.6;  1.0 + Thtv0];

%% Run Simulation

% clear storage variable
clear Wh_dat Omh_dat Thth_dat Wv_dat Omv_dat Thtv_dat uh_dat uv_dat ti

% Simulation Loop
for i = 1:Sim_samples
    % Store state vector values  
    Wh_dat(i) = x(1); Omh_dat(i) = x(2); Thth_dat(i) = x(3);
    Wv_dat(i) = x(4); Omv_dat(i) = x(5); Thtv_dat(i) = x(6);
    
    % Compute Reference
    TthtvRef = TthtvRef_v(i);
    TththRef = TththRef_v(i);
    [WhRef,OmhRef,WvRef,OmvRef] = compute_ref(TththRef,x(3),TthtvRef,x(6));
    ref = [WhRef OmhRef TththRef WvRef OmvRef TthtvRef-Thtv0]';
    
    % Adjust Vertical Angle State for MPC
    x_mpc = [x(1); x(2); x(3); x(4); x(5); x(6)-Thtv0];
    
    tic
    
    % --- LPV TRAJECTORY ESTIMATION METHODS ---
    % Choose only one method by uncommenting
    
    % Method 1: Frozen trajectory
    % Pk = compute_schedul_frozen(mpc, x_mpc, my_sched_fun);
    
    % Method 2: Iterative Fast trajectory (Warm-Start)
    Pk = compute_schedul_iterative_fast(mpc, x_mpc, my_sched_fun, n_rho);
    
    % Method 3: Iterative (SQP-like) trajectory refinement
    % Pk = compute_schedul_iterative(mpc, x_mpc, u_prev, ref, [], [], my_sched_fun, compute_A, compute_B, [], n_rho, n_iter);
    
    % Method 4: Recursive extrapolation trajectory
    % Pk = compute_schedul_recursive(mpc, x_mpc, n_rho, my_sched_fun, my_jacob_fun, rho_min, rho_max);
    
    % Build 3D affine matrix arrays using the interface
    A_lpv = traj_mat(compute_A, Pk, n_rho, n);
    B_lpv = traj_mat(compute_B, Pk, n_rho, n);
    
    % Update mpc problem structure with 3D arrays
    mpc = update_mpc_dynamics(mpc, A_lpv, B_lpv, []);
    
    % Solve mpc iteration
    [u_prev,iter,mpc] = mpc_solve(mpc, x_mpc, u_prev, ref, [], [], [], []);
    
    ti(i) = toc;
    
    % Store control actions
    uh_dat(i) = u_prev(1);
    uv_dat(i) = u_prev(2);
    
    % Run TRMS simulation & Forward Euler step
    dt_x = TRMS(x(1), x(2), x(3), x(4), x(5), x(6), u_prev(1), u_prev(2));
    x = x + Ts*dt_x;
end

%% Plots
figure

ax1 = subplot(2,2,1);
plot(time,TththRef_v)
hold on
plot(time,Thth_dat)
grid on
title('Horizontal Angle')
xlabel('Time (s)')
ylabel('Angle (rad)')
legend('Ref. \theta_h','\theta_h')
grid on

ax2 = subplot(2,2,2);
plot(time,TthtvRef_v-Thtv0)
hold on
plot(time,Thtv_dat-Thtv0)
title('Vertical Angle')
xlabel('Time (s)')
ylabel('Angle (rad)')
legend('Ref. \theta_v - \theta_{v0}','\theta_v - \theta_{v0}')
grid on

ax3 = subplot(2,2,3);
plot(time,uh_dat)
hold on
plot(time(1:end-1),diff(uh_dat))
grid on
title('Horizontal Fan Control Action')
xlabel('Time (s)')
ylabel('Fan Voltage (V)')
legend('u_h','\Delta u_h')
grid on

ax4 = subplot(2,2,4);
plot(time,uv_dat)
hold on
plot(time(1:end-1),diff(uv_dat))
grid on
title('Vertical Fan Control Action')
xlabel('Time (s)')
ylabel('Fan Voltage (V)')
legend('u_v','\Delta u_v')
grid on

linkaxes([ax1,ax3],'x')
linkaxes([ax2,ax4],'x')

figure
plot(time,ti)
title('Compute Time (s)')
xlabel('Time (s)')
grid on

%% Local functions for discrete matrices
function Ad = get_discrete_A(rho, Ts)
    % Extracts continuous matrix and discretizes it
    sys = qLPV_TRMS_SS(rho(1), rho(2), rho(3), rho(4), rho(5));
    [Ad, ~, ~] = init_discretize_system(sys.A, sys.B, [], Ts, 'forward');
end

function Bd = get_discrete_B(rho, Ts)
    % Extracts continuous matrix and discretizes it
    sys = qLPV_TRMS_SS(rho(1), rho(2), rho(3), rho(4), rho(5));
    [~, Bd, ~] = init_discretize_system(sys.A, sys.B, [], Ts, 'forward');
end