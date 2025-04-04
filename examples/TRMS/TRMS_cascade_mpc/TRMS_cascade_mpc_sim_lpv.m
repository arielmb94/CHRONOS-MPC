%% Call the mpc problem initialization script

TRMS_cascade_mpc_init

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

% Define masks to extract the individual controls actions sequences from
% the MIMO mpc
mask_Wh = 1:mpc.nu:mpc.Nu;
mask_Wv = 2:mpc.nu:mpc.Nu;

%% Run Simulation

clear Wh_dat Omh_dat Thth_dat Wv_dat Omv_dat Thtv_dat uh_dat uv_dat ti ...
    WhRef_dat WvRef_dat OmhRef_dat OmvRef_dat

% Simulation Loop
for i = 1:Sim_samples

% Assign state vector variables    
Wh   = x(1);    % Horizontal Fan Angular Speed
Omh  = x(2);    % Horizontal Angular Rate
Thth = x(3);    % Horizontal Angle
Wv   = x(4);    % Vertical Fan Angular Speed
Omv  = x(5);    % Vertical Angular Rate
Thtv = x(6);    % Vertical Angle

% Store state vector values for plotting and analysis  
Wh_dat(i)   = x(1);
Omh_dat(i)  = x(2);
Thth_dat(i) = x(3);
Wv_dat(i)   = x(4);
Omv_dat(i)  = x(5);
Thtv_dat(i) = x(6);

% Assign current reference values
TthtvRef = TthtvRef_v(i);
TththRef = TththRef_v(i);

% Compute Reference for rotors angular speed
OmhRef = (TththRef-Thth)/0.5;
OmvRef = (TthtvRef-Thtv)/0.5;

t0 = cputime;
% Update LPV model to current scheduling values
[A,B,Bd,Ah,Bh,Av,Bv] = qLPV_TRMS_cascade_mpc_SS(Wh,Omh,Thth,Wv,Thtv);
% Update MIMO mpc problem structure
mpc = update_mpc_sys_dynamics(mpc,eye(4)+Ts*A,Ts*B,Ts*Bd);
% Update Horizontal Fan mpc problem structure
mpc_h = update_mpc_sys_dynamics(mpc_h,1+Ts*Ah,Ts*Bh,[]);
% Update Vertical Fan mpc problem structure
mpc_v = update_mpc_sys_dynamics(mpc_v,1+Ts*Av,Ts*Bv,[]);

% Adjust Vertical Angle State
x_mpc = [Omh;Thth;Omv;Thtv-Thtv0];
% Define reference vector for MIMO mpc
ref = [OmhRef TththRef OmvRef TthtvRef-Thtv0]';

% solve MIMO mpc
[u_prev,J,x0] = mpc_solve(x0,x_mpc,u_prev,ref,uv,mpc,[],[],[]);

% Extract control action sequence from the optimization vector
W_ref = get_u(x0,mpc.nx,mpc.nu,mpc.N_ctr_hor,mpc.Nu);
% Use masks to isolate the Horizontal and Vertical fan speed references
% computed by the MIMO mpc
WhRef = W_ref(mask_Wh);
WvRef = W_ref(mask_Wv);

% solve Horizontal Fan mpc
[uh,J,x0_h] = mpc_solve(x0_h,Wh,uh,WhRef(1:mpc_h.Nu-1),[],mpc_h,WhRef(end),[],[]);
% solve Vertical Fan mpc
[uv,J,x0_v] = mpc_solve(x0_v,Wv,uv,WvRef(1:mpc_v.Nu-1),[],mpc_v,WvRef(end),[],[]);
ti(i) = cputime-t0;

% Storoge control action values for plotting and analysis
WhRef_dat(i) = WhRef(1);
WvRef_dat(i) = WvRef(1);
uh_dat(i) = uh;
uv_dat(i) = uv;

% Run TRMS simulation
dt_x = TRMS(Wh,Omh,Thth,Wv,Omv,Thtv,uh,uv);
% Forward euler step
x = x + Ts*dt_x;

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
plot(time,uh_dat)
hold on
plot(time(1:end-1),diff(uh_dat))
grid on
title('Horizontal Fan Control Action')
xlabel('Time (s)')
ylabel('Fan Voltage (V)')
legend('u_h','\Delta u_h')
grid on

ax4 = subplot(3,2,4);
plot(time,uv_dat)
hold on
plot(time(1:end-1),diff(uv_dat))
grid on
title('Vertical Fan Control Action')
xlabel('Time (s)')
ylabel('Fan Voltage (V)')
legend('u_v','\Delta u_v')
grid on

ax5 = subplot(3,2,5);
plot(time,WhRef_dat)
hold on
plot(time(1:end-1),diff(WhRef_dat))
grid on
title('MPC Computed Horizontal Fan Speed Reference')
xlabel('Time (s)')
ylabel('Fan Voltage (V)')
legend('\omega_h^{ref}','\Delta \omega_h^{ref}')
grid on

ax6 = subplot(3,2,6);
plot(time,WvRef_dat)
hold on
plot(time(1:end-1),diff(WvRef_dat))
grid on
title('MPC Computed Vertcal Fan Speed Reference')
xlabel('Time (s)')
ylabel('Fan Speed (V)')
legend('\omega_v^{ref}','\Delta \omega_v^{ref}')
grid on

linkaxes([ax1,ax3,ax5],'x')
linkaxes([ax2,ax4,ax6],'x')

figure
plot(time,ti)
title('Compute Time (s)')
xlabel('Time (s)')
grid on

