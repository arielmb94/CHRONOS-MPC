clear all
TRMS_init


TththRef = 0.0;
TthtvRef = 0.3;

%%
clear Wh_dat Omh_dat Thth_dat Wv_dat Omv_dat Thtv_dat uh_dat uv_dat

for i = 1:1000

Wh   = x(1);
Omh  = x(2);
Thth = x(3);
Wv   = x(4);
Omv  = x(5);
Thtv = x(6);

Wh_dat(i)   = x(1);
Omh_dat(i)  = x(2);
Thth_dat(i) = x(3);
Wv_dat(i)   = x(4);
Omv_dat(i)  = x(5);
Thtv_dat(i) = x(6);

OmhRef = (TththRef-Thth)/3;
OmvRef = (TthtvRef-(Thtv-Thtv0))/5;

ref = [Wh OmhRef TththRef Wv OmvRef TthtvRef]';

sys = qLPV_TRMS_SS(Wh,Omh,Thth,Wv,Thtv);
mpc = update_mpc_sys_dynamics(mpc,eye(6)+Ts*sys.A,Ts*sys.B,[]);

% Adjust state 6
x_mpc = [Wh;Omh;Thth;Wv;Omv;Thtv-Thtv0];


[u_prev,J,x0] = mpc_solve(x0,x_mpc,u_prev,ref,[],mpc,[],[],[]);


uh = u_prev(1);
uv = u_prev(2);

uh_dat(i) = u_prev(1);
uv_dat(i) = u_prev(2);


dt_x = TRMS(Wh,Omh,Thth,Wv,Omv,Thtv,uh,uv);
x = x + Ts*dt_x;

u_prev;

end

%%
close all
plot(Thth_dat)
grid on

figure
plot(Thtv_dat-Thtv0)
grid on

figure
plot(uh_dat)
hold on
plot(diff(uh_dat))
grid on

figure
plot(uv_dat)
hold on
plot(diff(uv_dat))
grid on

