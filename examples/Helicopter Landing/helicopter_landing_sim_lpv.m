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


t1 = cputime;
for i = 1:30
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
    
    xref = [r_x;r_vx;r_y;r_vy];
    
    [u_prev,J,x0] = mpc_solve(x0,x_prev,u_prev,xref,d,mpc,[],[],[]);
    
    px_dot = vx;
    vx_dot = u_prev(1) + kx * vx^2;
    py_dot = vy;
    vy_dot = u_prev(2) + ky * vy^2 + g;
    
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
end
t2 = cputime-t1;
t2/2400

%%
close all

figure
subplot(4,1,1)
plot(px_log)

subplot(4,1,2)
plot(vx_log)

subplot(4,1,3)
plot(py_log)

subplot(4,1,4)
plot(vy_log)

figure
subplot(2,1,1)
plot(ax_log)

subplot(2,1,2)
plot(ay_log)