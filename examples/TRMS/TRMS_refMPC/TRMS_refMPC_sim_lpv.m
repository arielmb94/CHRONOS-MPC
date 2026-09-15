%% Call the mpc problem initialization script

TRMS_refMPC_init

%% Define simulation duration and reference parameters

% Duration
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
TththRef_v = offset_TththRef + ampl_TththRef*sin(2*pi*(freq_TththRef)*time);
TthtvRef_v = offset_TthtvRef + ampl_TthtvRef*sin(2*pi*(freq_TthtvRef)*time);

% Initial values for fan speeds references
WhRef = x(1);
WvRef = x(4);

%% Run simulation

% Preallocate simulation data
Wh_dat = zeros(1,Sim_samples);
Omh_dat = zeros(1,Sim_samples);
Thth_dat = zeros(1,Sim_samples);
Wv_dat = zeros(1,Sim_samples);
Omv_dat = zeros(1,Sim_samples);
Thtv_dat = zeros(1,Sim_samples);
uh_dat = zeros(1,Sim_samples);
uv_dat = zeros(1,Sim_samples);
WhRef_dat = zeros(1,Sim_samples);
WvRef_dat = zeros(1,Sim_samples);
du_mpc_dat = zeros(4,Sim_samples);
t_dat = zeros(1,Sim_samples);
iter_dat = zeros(1,Sim_samples);

% Closed-loop simulation
for i = 1:Sim_samples

    % 1. Read the current plant state
    Wh = x(1);
    Omh = x(2);
    Thth = x(3);
    Wv = x(4);
    Omv = x(5);
    Thtv = x(6);

    % 2. Build the four-state body reference
    TththRef = TththRef_v(i);
    TthtvRef = TthtvRef_v(i);
    OmhRef = (TththRef-Thth)/0.5;
    OmvRef = (TthtvRef-Thtv)/0.5;
    ref = [OmhRef;TththRef;OmvRef;TthtvRef-Thtv0];

    % The terminal cost needs all six states. Reuse the virtual rotor
    % references optimized at the previous sample; no feedforward is needed.
    x_ref = [WhRef;OmhRef;TththRef;WvRef;OmvRef;TthtvRef-Thtv0];

    tic;

    % 3. Update the modified LPV model at the measured state
    sys = qLPV_TRMS_refMPC_SS(Wh,Omh,Thth,Wv,Thtv);
    mpc = update_mpc_dynamics(mpc,eye(6)+Ts*sys.A,Ts*sys.B,[]);

    % 4. Optimize two physical voltages and two virtual rotor references
    x_mpc = [Wh;Omh;Thth;Wv;Omv;Thtv-Thtv0];
    [u_mpc,mpc,iter] = mpc_solve(mpc,x_mpc,u_prev,ref,x_ref,[],[],[]);
    t_dat(i) = toc;

    uh = u_mpc(1);
    uv = u_mpc(2);
    WhRef = u_mpc(3);
    WvRef = u_mpc(4);

    % Store states, virtual references, and MPC actions
    Wh_dat(i) = Wh;
    Omh_dat(i) = Omh;
    Thth_dat(i) = Thth;
    Wv_dat(i) = Wv;
    Omv_dat(i) = Omv;
    Thtv_dat(i) = Thtv;
    uh_dat(i) = uh;
    uv_dat(i) = uv;
    WhRef_dat(i) = WhRef;
    WvRef_dat(i) = WvRef;
    du_mpc_dat(:,i) = u_mpc-u_prev;
    iter_dat(i) = iter;

    % 5. Apply only the physical voltage actions to the nonlinear plant
    dt_x = TRMS(Wh,Omh,Thth,Wv,Omv,Thtv,uh,uv);
    x = x + Ts*dt_x;
    u_prev = u_mpc;

end

%% Plots
figure

ax1 = subplot(3,2,1);
plot(time,TththRef_v)
hold on
plot(time,Thth_dat)
grid on
title('Horizontal Angle')
xlabel('Time (s)')
ylabel('Angle (rad)')
legend('Ref. \theta_h','\theta_h')
grid on

ax2 = subplot(3,2,2);
plot(time,TthtvRef_v-Thtv0)
hold on
plot(time,Thtv_dat-Thtv0)
title('Vertical Angle')
xlabel('Time (s)')
ylabel('Angle (rad)')
legend('Ref. \theta_v - \theta_{v0}','\theta_v - \theta_{v0}')
grid on

ax3 = subplot(3,2,3);
plot(time,WhRef_dat)
hold on
plot(time,Wh_dat)
grid on
title('Tail-Rotor Virtual Reference')
xlabel('Time (s)')
ylabel('Angular Speed (rad/s)')
legend('\omega_h^{ref}','\omega_h')
grid on

ax4 = subplot(3,2,4);
plot(time,WvRef_dat)
hold on
plot(time,Wv_dat)
grid on
title('Main-Rotor Virtual Reference')
xlabel('Time (s)')
ylabel('Angular Speed (rad/s)')
legend('\omega_v^{ref}','\omega_v')
grid on

ax5 = subplot(3,2,5);
plot(time,uh_dat,time,du_mpc_dat(1,:))
grid on
title('Tail-Rotor Control Action')
xlabel('Time (s)')
ylabel('Motor Voltage (V)')
legend('u_h','\Delta u_h')
grid on

ax6 = subplot(3,2,6);
plot(time,uv_dat,time,du_mpc_dat(2,:))
grid on
title('Main-Rotor Control Action')
xlabel('Time (s)')
ylabel('Motor Voltage (V)')
legend('u_v','\Delta u_v')
grid on

linkaxes([ax1,ax3,ax5],'x')
linkaxes([ax2,ax4,ax6],'x')

figure
subplot(2,1,1)
plot(time,t_dat)
ylabel('Online Controller Time (s)')
grid on

subplot(2,1,2)
stairs(time,iter_dat)
xlabel('Time (s)')
ylabel('Iterations')
grid on
