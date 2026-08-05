%% Call the mpc problem initialization script
stirring_tank_init
mpc.max_iter = 3;

%% Define simulation duration and reference parameters

% Duration
tsim = 240; % seconds
Sim_samples = tsim/Ts;
time = 0:Ts:tsim-Ts;

% Define concentration state reference
ref_c_vec = zeros(1,Sim_samples);
ref_c_vec(time < 90) = 0.27 + (0.65 - 0.27) * time(time < 90) / 90;
ref_c_vec(time >= 90 & time < 180) = 0.65;
ref_c_vec(time >= 180) = 0.65 - (0.65 - 0.3) * (time(time >= 180) - 180) / 60;

n = mpc.N;
n_rho = 3;
n_iter = 5;

% Setup scheduling functions and bounds
my_sched_fun = @(x) calc_rho_stirring_tank(x, k, M);
my_jacob_fun = @(x) calc_jacob_stirring_tank(x, k, M);

% LPV Matrix Builders
compute_A  = @(rho) eye(2) + Ts * [-1/theta_f - rho(1), -rho(2); 
                                    rho(1), -1/theta_f];
compute_B  = @(rho) Ts * [0; 
                          -alpha*(rho(3) - xc)];
compute_Bd = @(rho) Ts * [1/theta_f, rho(2); 
                          xf/theta_f, 0];

% Physical bounds for recursive method
x_min = [0.0; 0.1];
x_max = [1.0; 100.0];
rho_min = my_sched_fun(x_min);
rho_max = my_sched_fun(x_max);

%% Run Simulation

% Preallocate simulation data
c_dat = zeros(1,Sim_samples);
v_dat = zeros(1,Sim_samples);
u_dat = zeros(1,Sim_samples);
du_dat = zeros(1,Sim_samples);
ref_v_dat = zeros(1,Sim_samples);
t_dat = zeros(1,Sim_samples);
iter_dat = zeros(1,Sim_samples);

% Closed-loop simulation
for i = 1:Sim_samples
    % Update concentration state reference
    ref_c = ref_c_vec(i);
    
    % Compute temperature state reference from steady state equilibrium equation
    ref_v = - M/(log(1/(theta_f*k*ref_c)*(1-ref_c)));
    
    % Assign state vector variables  
    ck = x_prev(1);
    vk = x_prev(2);
    
    tic;
    
    % Update mpc disturbance vector
    d = [1;vk];
    
    % Update reference vector
    xref = [ref_c;ref_v];
    
    % --- LPV TRAJECTORY ESTIMATION METHODS ---
    % Choose only one method to use by uncommenting it and commenting the others
    
    % Method 1: Frozen trajectory
    %Pk = compute_schedul_frozen(mpc, x_prev, my_sched_fun);
    
    % Method 2: Iterative Fast trajectory
    %Pk = compute_schedul_iterative_fast(mpc, x_prev, my_sched_fun, n_rho);
    
    % Method 3: Iterative (SQP-like) trajectory refinement
     Pk = compute_schedul_iterative(mpc, x_prev, u_prev, xref, [], d, my_sched_fun, ...
                                    compute_A, compute_B, compute_Bd, n_rho, n_iter);
    
    % Method 4: Recursive extrapolation trajectory
    %Pk = compute_schedul_recursive(mpc, x_prev, n_rho, my_sched_fun, my_jacob_fun, rho_min, rho_max);
    
    % Build 3D affine matrix arrays using the interface
    A_lpv  = traj_mat(compute_A, Pk, n_rho, n);
    B_lpv  = traj_mat(compute_B, Pk, n_rho, n);
    Bd_lpv = traj_mat(compute_Bd, Pk, n_rho, n);
    
    % Update mpc problem dynamics
    mpc = update_mpc_dynamics(mpc, A_lpv, B_lpv, Bd_lpv);
    
    % Solve mpc iteration
    [u_prev, iter, mpc] = mpc_solve(mpc, x_prev, u_prev, xref, [], d, [], []);
    
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

figure
ax1 = subplot(3,1,1);
plot(time,ref_c_vec,'r',time,c_dat,'b')
grid on
ylim([0.2 0.7])
xlim([0 tsim])
legend('Concentration Ref.','c_k')
xlabel('Time (s)')
ylabel('Tank Concentration')

ax2 = subplot(3,1,2);
plot(time,ref_v_dat,'r',time,v_dat,'b')
grid on
ylim([0.5 0.7])
xlim([0 tsim])
legend('Temperature Ref.','v_k')
xlabel('Time (s)')
ylabel('Tank Temperature')

ax3 = subplot(3,1,3);
plot(time,u_dat,time,du_dat)
grid on
xlim([0 tsim])
legend('u','\Delta u')
xlabel('Time (s)')
ylabel('Coolant Flow Rate')

linkaxes([ax1,ax2,ax3],'x')

figure
subplot(2,1,1)
plot(time,t_dat)
xlim([0 tsim])
ylabel('Online Controller Time (s)')
grid on

subplot(2,1,2)
stairs(time,iter_dat)
xlim([0 tsim])
xlabel('Time (s)')
grid on

%% Local Functions
function rho = calc_rho_stirring_tank(x, k_cte, M)
    ck = x(1);
    vk = max(x(2), 0.1); % Avoid division by zero
    
    rho1 = k_cte * exp(-M/vk);
    rho2 = k_cte * ck * M * exp(-M/vk) / (vk^2);
    rho3 = vk;
    
    rho = [rho1; rho2; rho3];
end

function sigma_k = calc_jacob_stirring_tank(x, k_cte, M)
    ck = x(1);
    vk = max(x(2), 0.1); % Avoid division by zero
    
    d_rho1_dc = 0;
    d_rho1_dv = k_cte * M / (vk^2) * exp(-M/vk);
    
    d_rho2_dc = k_cte * M / (vk^2) * exp(-M/vk);
    d_rho2_dv = k_cte * ck * M * exp(-M/vk) * (M/(vk^4) - 2/(vk^3));
    
    d_rho3_dc = 0;
    d_rho3_dv = 1;
    
    sigma_k = [d_rho1_dc, d_rho1_dv;
               d_rho2_dc, d_rho2_dv;
               d_rho3_dc, d_rho3_dv];
end
