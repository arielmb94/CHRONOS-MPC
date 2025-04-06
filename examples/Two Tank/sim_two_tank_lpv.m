%% Call the mpc problem initialization script

two_tank_init

%% Define simulation duration and reference parameters

% Duration
tsim = 5; % seconds
Sim_samples = tsim/Ts;
time = 0:Ts:tsim-Ts;

% Define step reference
clear r
r = zeros(1,Sim_samples);
r(time<=2.5) = 0.7;
r(time>2.5) = 0.25;

% To avoid feasibility problems due to large step changes it
% is better to low-pass step references
tau = 0.1;      % time constant for reference filter
xf = h2;        % initial value for reference filter state

%% Run Simulation

% clear storage variables
clear rf_dat rk h1_dat h2_dat u_dat t_dat

% Simulation Loop
for k = 1:Sim_samples

% Assign state vector variables    
h1  = x_prev(1);        % Tank 1 water height
h2  = x_prev(2);        % Tank 2 water height

% Low pass reference filter step
xf = xf + Ts*(-xf/tau+r(k)/tau);
% Tracking vector for terminal constraint
x_ref = [xf;xf];

tic;
% Update LPV model
A_lpv = eye(2)+Ts*[-sqrt(2*g)*sqrt(h1)/(Ab*h1) 0;
     sqrt(2*g)*sqrt(h1)/(Ab*h1) -sqrt(2*g)*sqrt(h2)/(Ab*h2)];
% Update mpc problem dynamics
mpc = update_mpc_sys_dynamics(mpc,A_lpv,mpc.B,[]);

% Solve mpc iteration
[u_prev,J,x0] = mpc_solve(x0,x_prev,u_prev,xf,[],mpc,x_ref,[],[]);
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
ylabel('Water Height')

ax2 = subplot(2,1,2);
plot(time,u_dat,time(1:end-1),diff(u_dat))
grid on
legend('u','\Delta u')
xlabel('Time (s)')
ylabel('Input Mass Flow')

linkaxes([ax1,ax2 ],'x')

figure
plot(time,t_dat)
title('Compute Time (s)')
xlabel('Time (s)')
grid on
