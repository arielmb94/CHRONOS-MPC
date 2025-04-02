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

%% Run Simulation

% clear storage variable
clear Wh_dat Omh_dat Thth_dat Wv_dat Omv_dat Thtv_dat uh_dat uv_dat ti

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

% Compute Reference for each state
[WhRef,OmhRef,WvRef,OmvRef] = compute_ref(TththRef,Thth,TthtvRef,Thtv);

% Define reference vector
ref = [WhRef OmhRef TththRef WvRef OmvRef TthtvRef-Thtv0]';

t0 = cputime;
% Update LPV model to current scheduling values
sys = qLPV_TRMS_SS(Wh,Omh,Thth,Wv,Thtv);
% Update mpc problem structure
mpc = update_mpc_sys_dynamics(mpc,eye(6)+Ts*sys.A,Ts*sys.B,[]);

% Adjust Vertical Angle State
x_mpc = [Wh;Omh;Thth;Wv;Omv;Thtv-Thtv0];
% Solve mpc iteration
[u_prev,J,x0] = mpc_solve(x0,x_mpc,u_prev,ref,[],mpc,[],[],[]);
ti(i) = cputime-t0;

% Assign control actions
uh = u_prev(1);
uv = u_prev(2);
% Storoge control action values for plotting and analysis
uh_dat(i) = u_prev(1);
uv_dat(i) = u_prev(2);

% Run TRMS simulation
dt_x = TRMS(Wh,Omh,Thth,Wv,Omv,Thtv,uh,uv);
% Forward euler step
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
