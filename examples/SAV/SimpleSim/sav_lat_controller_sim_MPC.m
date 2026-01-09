%% Call the mpc problem initialization script

sav_lat_controller_init

%% Define simulation duration and reference parameters

% Duration
tsim = 10; % seconds
Sim_samples = tsim/Ts;
time = 0:Ts:tsim-Ts;

% Define step reference
clear r
r = zeros(1,Sim_samples);
r(time<=5) = 0.4;
r(time>5) = -0.4;

% To avoid feasibility problems due to large step changes it
% is better to low-pass step references
tau = 0.1;      % time constant for reference filter
xf = 0;        % initial value for reference filter state

delta = 0;
%% Run Simulation

% clear storage variables
clear rf_dat rk h1_dat h2_dat u_dat t_dat

% Simulation Loop
for k = 1:Sim_samples

% Assign state vector variables    
vy  = x_prev(1);        
yaw_rate  = x_prev(2);        

% Low pass reference filter step
xf = xf + Ts*(-xf/tau+r(k)/tau);

tic;

A = update_BM(vx);
Ak = eye(2)+Ts*A;
mpc = update_mpc_sys_dynamics(mpc,Ak,[],[]);

[delta,x0] = mpc_solve(mpc,x0,x_prev,delta,xf,[],[],[],[]);
tk = toc;

% Store variables values for plotting and analysis  
rf_dat(:,k) = xf;
vy_dat(:,k) = vy;
yaw_rate_dat(:,k) = yaw_rate;
u_dat(k) = delta;
t_dat(k) = tk;

% Forward Euler step of Two Tank nonlinear dynamics
x_dot = simNLPVCalpha(0,0,delta,[vx;x_prev],x);
x_dot = x_dot(2:3);

x_prev = x_prev+Ts*x_dot;

end

%% Plots

figure
ax1 = subplot(2,1,1);
plot(time,r,'r',time,rf_dat,'--r',time,yaw_rate_dat,'g')
grid on
legend('Reference','Filtered Reference','yaw rate')
xlabel('Time (s)')
ylabel('Water Height')