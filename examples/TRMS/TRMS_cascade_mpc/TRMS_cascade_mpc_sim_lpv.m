%% Call the MPC problem initialization script

TRMS_cascade_mpc_init

%% Define simulation duration and reference parameters

tsim = 200; % seconds
Sim_samples = tsim/Ts;
time = 0:Ts:tsim-Ts;

% Horizontal-angle sinusoidal reference parameters
freq_TththRef = 1/31; % Hz
offset_TththRef = 0;
ampl_TththRef = 1;

% Vertical-angle sinusoidal reference parameters
freq_TthtvRef = 1/47; % Hz
offset_TthtvRef = -0.6;
ampl_TthtvRef = 0.5;

% Sinusoidal angle references
TththRef_v = offset_TththRef + ampl_TththRef*sin(2*pi*freq_TththRef*time);
TthtvRef_v = offset_TthtvRef + ampl_TthtvRef*sin(2*pi*freq_TthtvRef*time);

%% Run simulation

% Preallocate simulation data
Wh_dat = zeros(1,Sim_samples);
Omh_dat = zeros(1,Sim_samples);
Thth_dat = zeros(1,Sim_samples);
Wv_dat = zeros(1,Sim_samples);
Omv_dat = zeros(1,Sim_samples);
Thtv_dat = zeros(1,Sim_samples);
WhRef_dat = zeros(1,Sim_samples);
WvRef_dat = zeros(1,Sim_samples);
uh_dat = zeros(1,Sim_samples);
uv_dat = zeros(1,Sim_samples);
duh_dat = zeros(1,Sim_samples);
duv_dat = zeros(1,Sim_samples);
t_dat = zeros(1,Sim_samples);
iter_outer_dat = zeros(1,Sim_samples);
iter_h_dat = zeros(1,Sim_samples);
iter_v_dat = zeros(1,Sim_samples);

% Closed-loop simulation
for i = 1:Sim_samples

    % 1. Read the current plant state
    Wh = x(1);
    Omh = x(2);
    Thth = x(3);
    Wv = x(4);
    Omv = x(5);
    Thtv = x(6);

    % 2. Build the outer-MPC reference from the requested angles
    TththRef = TththRef_v(i);
    TthtvRef = TthtvRef_v(i);
    OmhRef = (TththRef-Thth)/0.5;
    OmvRef = (TthtvRef-Thtv)/0.5;

    x_outer = [Omh;Thth;Omv;Thtv-Thtv0];
    ref_outer = [OmhRef;TththRef;OmvRef;TthtvRef-Thtv0];

    tic

    % 3. Update the outer and inner LPV models at the measured state
    [A,B,Bd,Ah,Bh,Av,Bv] = qLPV_TRMS_cascade_mpc_SS(Wh,Omh,Thth,Wv,Thtv);
    mpc = update_mpc_dynamics(mpc,eye(4)+Ts*A,Ts*B,Ts*Bd);
    mpc_h = update_mpc_dynamics(mpc_h,1+Ts*Ah,Ts*Bh,[]);
    mpc_v = update_mpc_dynamics(mpc_v,1+Ts*Av,Ts*Bv,[]);

    % 4. The outer MPC computes both rotor-reference sequences. The current
    % main-rotor voltage enters its model through the disturbance input.
    [omega_ref_k,mpc,iter_outer] = mpc_solve(mpc,x_outer,omega_ref_prev,ref_outer,[],uv_prev,[],[]);
    WhRef_seq = mpc.u(1,:);
    WvRef_seq = mpc.u(2,:);

    % 5. Each inner MPC tracks the complete sequence from the outer MPC
    [uh_k,mpc_h,iter_h] = mpc_solve(mpc_h,Wh,uh_prev,WhRef_seq,WhRef_seq(end),[],[],[]);
    [uv_k,mpc_v,iter_v] = mpc_solve(mpc_v,Wv,uv_prev,WvRef_seq,WvRef_seq(end),[],[],[]);
    t_dat(i) = toc;

    % Store states, first rotor references, and applied voltages
    Wh_dat(i) = Wh;
    Omh_dat(i) = Omh;
    Thth_dat(i) = Thth;
    Wv_dat(i) = Wv;
    Omv_dat(i) = Omv;
    Thtv_dat(i) = Thtv;
    WhRef_dat(i) = WhRef_seq(1);
    WvRef_dat(i) = WvRef_seq(1);
    uh_dat(i) = uh_k;
    uv_dat(i) = uv_k;
    duh_dat(i) = uh_k-uh_prev;
    duv_dat(i) = uv_k-uv_prev;
    iter_outer_dat(i) = iter_outer;
    iter_h_dat(i) = iter_h;
    iter_v_dat(i) = iter_v;

    % 6. Apply the inner MPC voltages to the nonlinear plant
    dt_x = TRMS(Wh,Omh,Thth,Wv,Omv,Thtv,uh_k,uv_k);
    x = x + Ts*dt_x;

    omega_ref_prev = omega_ref_k;
    uh_prev = uh_k;
    uv_prev = uv_k;

end

%% Plots

figure

ax1 = subplot(3,2,1);
plot(time,TththRef_v,time,Thth_dat)
grid on
title('Horizontal Angle')
xlabel('Time (s)')
ylabel('Angle (rad)')
legend('Ref. \theta_h','\theta_h')

ax2 = subplot(3,2,2);
plot(time,TthtvRef_v-Thtv0,time,Thtv_dat-Thtv0)
grid on
title('Vertical Angle')
xlabel('Time (s)')
ylabel('Angle (rad)')
legend('Ref. \theta_v - \theta_{v0}','\theta_v - \theta_{v0}')

ax3 = subplot(3,2,3);
plot(time,WhRef_dat,time,Wh_dat)
grid on
title('Tail-Rotor Speed')
xlabel('Time (s)')
ylabel('Angular Speed (rad/s)')
legend('\omega_h^{ref}','\omega_h')

ax4 = subplot(3,2,4);
plot(time,WvRef_dat,time,Wv_dat)
grid on
title('Main-Rotor Speed')
xlabel('Time (s)')
ylabel('Angular Speed (rad/s)')
legend('\omega_v^{ref}','\omega_v')

ax5 = subplot(3,2,5);
plot(time,uh_dat,time,duh_dat)
grid on
title('Tail-Rotor Control Action')
xlabel('Time (s)')
ylabel('Motor Voltage (V)')
legend('u_h','\Delta u_h')

ax6 = subplot(3,2,6);
plot(time,uv_dat,time,duv_dat)
grid on
title('Main-Rotor Control Action')
xlabel('Time (s)')
ylabel('Motor Voltage (V)')
legend('u_v','\Delta u_v')

linkaxes([ax1,ax3,ax5],'x')
linkaxes([ax2,ax4,ax6],'x')

figure
subplot(2,1,1)
plot(time,t_dat)
ylabel('Online Controller Time (s)')
grid on

subplot(2,1,2)
stairs(time,iter_outer_dat)
hold on
stairs(time,iter_h_dat)
stairs(time,iter_v_dat)
xlabel('Time (s)')
ylabel('Iterations')
legend('Outer MIMO','Tail rotor','Main rotor')
grid on
