% Load init script defining the mpc problem
helicopter_landing_init

clear c_dat v_dat uk rk
px = x_prev(1);
vx = x_prev(2);
py = x_prev(3);
vy = x_prev(4);

r_x = 300;
r_vx = 0;
r_y = 10;
r_vy = 0;

f_x = px;
f_vx = vx;
f_y = py;
f_vy = vy;
tau = 15;


t1 = cputime;
for i = 1:700
    f_x = f_x + Ts*(-f_x/tau+r_x/tau);
    f_vx = f_vx + Ts*(-f_vx/tau+r_vx/tau);
    f_y = f_y + Ts*(-f_y/tau+r_y/tau);
    f_vy = f_vy + Ts*(-f_vy/tau+r_vy/tau);

    % Update LPV model
    rho_1 = -kx * vx;
    rho_2 = -ky * vy;
    
    A_lpv = [1 1*Ts 0 0; 
            0 1+Ts*rho_1 0 0;
            0 0 1 1*Ts;
            0 0 0 1+Ts*rho_2];
    
    B_lpv = [0 0;
             Ts 0;
             0 0;
             0 Ts];
    
    Bd_lpv = [0;
        0;
        0;
        Ts];
    
    mpc = update_mpc_sys_dynamics(mpc,A_lpv,B_lpv,Bd_lpv);
    
    d = [g];
    
    xref = [f_x;f_vx;f_y;f_vy];
    
    [u_prev,J,x0] = mpc_solve(x0,x_prev,u_prev,xref,d,mpc,[],[],[]);
    
    px_dot = vx;
    vx_dot = u_prev(1) - kx * vx^2;
    py_dot = vy;
    vy_dot = u_prev(2) - ky * vy^2 + g;
    
    px = px + Ts*px_dot;
    vx = vx + Ts*vx_dot;
    py = py + Ts*py_dot;
    vy = vy + Ts*vy_dot;
    
    x_prev = [px;vx;py;vy];
    
    px_log(:,i) = px;
    py_log(:,i) = py;
    vx_log(:,i) = vx;
    vy_log(:,i) = vy;
    ax_log(i) = u_prev(1);
    ay_log(i) = u_prev(2);

    ref_x_log(:,i) = xref(1);
    ref_vx_log(:,i) = xref(2);
    ref_y_log(:,i) = xref(3);
    ref_vy_log(:,i) = xref(4);
end
t2 = cputime-t1;
t2/2400

%%
close all

figure
subplot(4,1,1)
plot(px_log)
hold on
plot(ref_x_log)
legend('px', 'ref')
ylabel("px")

subplot(4,1,2)
plot(vx_log)
hold on
plot(ref_vx_log)
legend('vx', 'ref')
ylabel("vx")

subplot(4,1,3)
plot(py_log)
hold on
plot(ref_y_log)
legend('py', 'ref')
ylabel("py")

subplot(4,1,4)
plot(vy_log)
hold on
plot(ref_vy_log)
legend('vy', 'ref')
ylabel("vy")

figure
subplot(2,1,1)
plot(ax_log)
ylabel("ax")

subplot(2,1,2)
plot(ay_log)
ylabel("ay")

figure
plot(px_log, py_log)
xlabel('X Axis')
ylabel('Y Axis')
title("2D Helicopter View")