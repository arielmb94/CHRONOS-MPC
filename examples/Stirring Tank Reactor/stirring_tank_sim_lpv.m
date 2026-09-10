%% Call the mpc problem initialization script
stirring_tank_init
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

%% Run simulation

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

    % 1. Read the current state and update the full-state reference
    ck = x_prev(1);
    vk = x_prev(2);
    ref_c = ref_c_vec(i);
    ref_v = -M/log((1-ref_c)/(theta_f*k*ref_c));
    x_ref = [ref_c;ref_v];

    tic

    % 2. Evaluate and discretize the hybrid model at the measured state
    A_lpv = eye(2)+Ts*[-1/theta_f-k*exp(-M/vk) -k*ck*M*exp(-M/vk)/(vk^2);
         k*exp(-M/vk) -1/theta_f];
    B_lpv = Ts*[0; -alpha*(vk-xc)];
    Bd_lpv = Ts*[1/theta_f k*ck*M*exp(-M/vk)/(vk^2); xf/theta_f 0];

    mpc = update_mpc_dynamics(mpc,A_lpv,B_lpv,Bd_lpv);
    d = [1;vk];  % constant and Taylor-expansion offset

    % 3. Solve with the tracking reference and known input. Terminal and
    % custom-signal inputs are not used in this example.
    [u_k,iter,mpc] = mpc_solve(mpc,x_prev,u_prev,x_ref,[],d,[],[]);
    tk = toc;

    % Keep the returned mpc for the next control sample
    c_dat(i) = ck;
    v_dat(i) = vk;
    u_dat(i) = u_k;
    du_dat(i) = u_k-u_prev;
    ref_v_dat(i) = ref_v;
    t_dat(i) = tk;
    iter_dat(i) = iter;

    % 4. Apply the first control action to the nonlinear plant
    c_next = ck + Ts*((1-ck)/theta_f - k*ck*exp(-M/vk));
    v_next = vk + Ts*((xf-vk)/theta_f + k*ck*exp(-M/vk)...
                      -alpha*u_k*(vk-xc));

    x_prev = [c_next;v_next];
    u_prev = u_k;

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
ylabel('Iterations')
grid on
