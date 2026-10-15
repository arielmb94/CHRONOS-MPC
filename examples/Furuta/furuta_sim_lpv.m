furuta_init

%% simulation
tsim = 25; % seconds
Sim_samples = tsim / Ts;
time = 0:Ts:tsim-Ts;

% iterative fast method
use_iterative_fast = true;

%% Preallocate simulation data
th1_dat  = zeros(1, Sim_samples);
th2_dat  = zeros(1, Sim_samples);
dth1_dat = zeros(1, Sim_samples);
dth2_dat = zeros(1, Sim_samples);
ref1_dat = zeros(1, Sim_samples);
u_dat    = zeros(1, Sim_samples);
du_dat   = zeros(1, Sim_samples);
t_dat    = zeros(1, Sim_samples);
iter_dat = zeros(1, Sim_samples);

Ac_lpv = zeros(4, 4, N);
Mc_lpv = zeros(4, 4, N);

%% Closed-loop simulation
for i = 1:Sim_samples
    i
    t_curr = time(i);
    
    % 1. Reference generation (Step on arm angle theta1 at t = 11 s)
    if t_curr >= 11.0
        th1_ref = pi/2; % Reference step
    else
        th1_ref = 0.0;
    end
    ref = [th1_ref; 0; 0; 0];
    
    % Input disturbance at t = 20 s (impulse/short pulse)
    if t_curr >= 20.0 && t_curr < 20.1
        d_u = 0.3;
    else
        d_u = 0.0;
    end
    
    tic
    % Implement iterative fast
    if use_iterative_fast && i > 1
        X_pred = mpc.s;
        
        % use current measured state
        [Ac_lpv(:,:,1), Bc, Mc_lpv(:,:,1)] = qLPV_Furuta_SS(x);
        for k = 2:N
            % use shifted predictions
            [Ac_lpv(:,:,k), ~, Mc_lpv(:,:,k)] = qLPV_Furuta_SS(X_pred(:, k-1));
        end
        mpc = update_mpc_dyn_continuous_forward(mpc, Ac_lpv, Bc, [], Ts, Mc_lpv);
    else
        % frozen
        [Ac, Bc, Mc] = qLPV_Furuta_SS(x);
        mpc = update_mpc_dyn_continuous_forward(mpc, Ac, Bc, [], Ts, Mc);
    end
    
    [u_k, mpc, iter] = mpc_solve(mpc, x, u_prev, ref, [], [], [], []);
    t_dat(i) = toc;
    
    % Store data
    th1_dat(i)  = x(1);
    th2_dat(i)  = x(2);
    dth1_dat(i) = x(3);
    dth2_dat(i) = x(4);
    ref1_dat(i) = th1_ref;
    u_dat(i)    = u_k;
    du_dat(i)   = u_k - u_prev;
    iter_dat(i) = iter;
    
    % nonlinear plant (RK4 for high-fidelity plant simulation)
    u_applied = u_k + d_u;
    k1 = Furuta(x, u_applied);
    k2 = Furuta(x + 0.5*Ts*k1, u_applied);
    k3 = Furuta(x + 0.5*Ts*k2, u_applied);
    k4 = Furuta(x + Ts*k3, u_applied);
    x = x + (Ts/6)*(k1 + 2*k2 + 2*k3 + k4);
    
    u_prev = u_k;
end

%% Plots
figure
ax1 = subplot(3,1,1);
plot(time, ref1_dat, '--k', 'LineWidth', 1.2); hold on;
plot(time, th1_dat, 'b', 'LineWidth', 1.2);
grid on; ylabel('\theta_1 (rad)');
title('(a) Arm Angle \theta_1');
legend('Reference', '\theta_1');

ax2 = subplot(3,1,2);
yline(0, '--k', 'LineWidth', 1.2); hold on;
plot(time, th2_dat, 'r', 'LineWidth', 1.2);
grid on; ylabel('\theta_2 (rad)');
title('(b) Pendulum Angle \theta_2 (0 = Upright)');

ax3 = subplot(3,1,3);
plot(time, u_dat, 'g', 'LineWidth', 1.2);
grid on; ylabel('u (N.m)'); xlabel('Time (s)');
title('(c) Control Input u');
ylim([-0.7 0.7]);
linkaxes([ax1, ax2, ax3], 'x');

figure
subplot(2,1,1)
semilogy(time, t_dat * 1e3, 'LineWidth', 1.2); hold on;
yline(Ts * 1e3, '--r', 'Ts (20 ms)');
ylabel('Compute Time (ms)');
title('Online Controller Computation Time');
grid on;

subplot(2,1,2)
stairs(time, iter_dat, 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('Iterations');
grid on;