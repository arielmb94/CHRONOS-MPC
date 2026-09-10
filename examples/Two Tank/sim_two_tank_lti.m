%% Call the mpc problem initialization script

two_tank_init

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

    % 1. Read the current state and update the references
    h1 = x_prev(1);
    h2 = x_prev(2);
    xf = xf + Ts*(-xf/tau+r(k)/tau);
    x_ref = [xf;xf];  % terminal target: equal steady-state tank heights

    tic;

    % 2. Solve using the current state, previous input, output reference, and
    % terminal-state reference. Empty arguments mean no d, dz, or dh inputs.
    [u_k,iter,mpc] = mpc_solve(mpc,x_prev,u_prev,xf,x_ref,[],[],[]);
    tk = toc;

    % Keep the returned mpc for the next control sample
    rf_dat(k) = xf;
    h1_dat(k) = h1;
    h2_dat(k) = h2;
    u_dat(k) = u_k;
    du_dat(k) = u_k-u_prev;
    t_dat(k) = tk;
    iter_dat(k) = iter;

    % 3. Apply the first control action to the nonlinear plant
    h1_next = h1 + Ts*(u_k-sqrt(2*g*h1))/Ab;
    h2_next = h2 + Ts*(sqrt(2*g*h1)-sqrt(2*g*h2))/Ab;

    x_prev = [h1_next;h2_next];
    u_prev = u_k;

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
