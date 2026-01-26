function y = sav_lat_mpc_solve(mpc,x0,vx,vy,yaw_rate,yaw_rate_ref,...
                                                    u_prev)

tic
[A,B] = update_BM(vx);
Ak = eye(2)+0.0200*A;
Bk = 0.02*B;
mpc = update_mpc_sys_dynamics(mpc,Ak,Bk,[]);

x_prev = [vy;yaw_rate];

[u,x0] = mpc_solve(mpc,x0,x_prev,u_prev,yaw_rate_ref,[],[],[],[]);

t = toc;

y = [u;t;x0];
end