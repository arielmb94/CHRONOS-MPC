%% Call the mpc problem initialization script

two_tank_init
mpc.max_iter = 3;
%% Define simulation duration and reference parameters

% Duration
tsim = 5; % seconds
Sim_samples = tsim/Ts;
time = 0:Ts:tsim-Ts;

% Define step reference
r = zeros(1,Sim_samples);
r(time<=2.5) = 0.7;
r(time>2.5) = 0.25;

% Smooth the commanded trajectory so that the controller does not react to
% the full reference step at once. This is useful with input-rate limits and
% a fixed, small Newton-iteration budget; it is not a feasibility mechanism.
tau = 0.1;      % time constant for reference filter
xf = h2;        % initial value for reference filter state

n = mpc.N;
n_rho = 2;
n_iter = 5;

% Setup scheduling functions and bounds
my_sched_fun = @(x) [1/sqrt(max(x(1), 1e-4)); 
                     1/sqrt(max(x(2), 1e-4))];
                     
my_jacob_fun = @(x) [-0.5/(max(x(1), 1e-4)^1.5), 0; 
                      0, -0.5/(max(x(2), 1e-4)^1.5)];

k_tank = sqrt(2*g)/Ab;
compute_Ac = @(rho) [-k_tank*rho(1), 0; 
                      k_tank*rho(1), -k_tank*rho(2)];
Bc = [1/Ab; 0];

rho_min = my_sched_fun([1.0; 1.0]); 
rho_max = my_sched_fun([0.01; 0.01]); 

A_lpv = zeros(2,2,n);

%% Run Simulation

% initialize varying A matrix
A_lpv = zeros(2,2,n);

%% Run simulation

% Preallocate simulation data
rf_dat = zeros(1,Sim_samples);
h1_dat = zeros(1,Sim_samples);
h2_dat = zeros(1,Sim_samples);
u_dat = zeros(1,Sim_samples);
du_dat = zeros(1,Sim_samples);
t_dat = zeros(1,Sim_samples);
iter_dat = zeros(1,Sim_samples);

% Closed-loop simulation
for k = 1:Sim_samples
    % Assign state vector variables    
    h1  = x_prev(1);        % Tank 1 water height
    h2  = x_prev(2);        % Tank 2 water height
    
    % Low pass reference filter step
    xf = xf + Ts*(-xf/tau+r(k)/tau);
    
    % Tracking vector for terminal constraint
    x_ref = [xf;xf];
    
    tic;
    
    % --- LPV TRAJECTORY ESTIMATION METHODS ---
    % Choose only one method to use by uncommenting it and commenting the others
    
    % Method 1: Frozen trajectory
    %Pk = compute_schedul_frozen(mpc, x_prev, my_sched_fun);
    
    % Method 2: Iterative Fast trajectory
    Pk = compute_schedul_iterative_fast(mpc, x_prev, my_sched_fun, n_rho);
    
    % Method 3: Iterative (SQP-like) trajectory refinement
    %Pk = compute_schedul_iterative(mpc, x_prev, u_prev, xf, x_ref, [], my_sched_fun, compute_A, [], [], n_rho, n_iter);
    
    % Method 4: Recursive extrapolation trajectory
    % Pk = compute_schedul_recursive(mpc, x_prev, n_rho, my_sched_fun, my_jacob_fun, rho_min, rho_max);
    
    % Build 3D affine matrix array for Ac using the interface
    Ac_lpv = traj_mat(compute_Ac, Pk, n_rho, n);
    
    % Update mpc problem dynamics using the continuous API
    mpc = update_mpc_dyn_continuous_forward(mpc, Ac_lpv, Bc, [], Ts);
    
    % Solve mpc iteration
    [u_prev, mpc, iter] = mpc_solve(mpc, x_prev, u_prev, xf, x_ref, [], [], []);
    
    tk = toc;
    
    % Store variables values for plotting and analysis  
    rf_dat(:,k) = xf;
    h1_dat(:,k) = h1;
    h2_dat(:,k) = h2;
    u_dat(k) = u_prev;
    t_dat(k) = tk;
    
    % Forward Euler step of Two Tank nonlinear dynamics
    h1 = h1 + Ts*(u_prev/Ab-sqrt(2*g)*sqrt(h1)/Ab);
    h2 = h2 + Ts*(sqrt(2*g)*sqrt(h1)/Ab-sqrt(2*g)*sqrt(h2)/Ab);
    
    % update state vector for the following iteration
    x_prev = [h1;h2];
end

%% Plots

figure
ax1 = subplot(2,1,1);
plot(time,r,'r',time,rf_dat,'--r',time,h1_dat,'g',time,h2_dat,'b')
grid on
legend('Reference','Filtered Reference','h1','h2')
xlabel('Time (s)')
ylabel('Water Height (m)')

ax2 = subplot(2,1,2);
plot(time,u_dat,time,du_dat)
grid on
legend('u','\Delta u')
xlabel('Time (s)')
ylabel('Inlet Flow (m^3/s)')

linkaxes([ax1,ax2 ],'x')

figure
subplot(2,1,1)
plot(time,t_dat)
ylabel('Compute Time (s)')
grid on

subplot(2,1,2)
stairs(time,iter_dat)
xlabel('Time (s)')
ylabel('Iterations')
grid on
