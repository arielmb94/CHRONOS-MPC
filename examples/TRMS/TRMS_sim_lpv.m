Thtv0 = -0.619426178368110;
x = [0.01,0.010,0.01,0.01,0.01,Thtv0+0.01];

Wh   = x(1);
Omh  = x(2);
Thth = x(3);
Wv   = x(4);
Omv = x(5); 
Thtv = x(6);
%%
TRMS_init
%%

tic
[u_prev,J,x0] = mpc_solve(x0,x,u_prev,x+[0 0 0.1 0 0 0.1]',[],mpc,[],[],[]);
toc

uh = u_prev(1)
uv = u_prev(2)


dt_x = TRMS(Wh,Omh,Thth,Wv,Omv,Thtv,uh,uv)
x = x+Ts*dt_x;

Wh   = x(1)
Omh  = x(2)
Thth = x(3)
Wv   = x(4)
Omv = x(5) 
Thtv = x(6)