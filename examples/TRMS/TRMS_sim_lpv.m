clear all
TRMS_init


TththRef = 0.3;
TthtvRef = -0.4;

%%
clear Wh_dat Omh_dat Thth_dat Wv_dat Omv_dat Thtv_dat uh_dat uv_dat ti

Sim_samples = 30/Ts;

for i = 1:Sim_samples

Wh   = x(1);
Omh  = x(2);
Thth = x(3);
Wv   = x(4);
Omv  = x(5);
Thtv = x(6);

[WhRef,OmhRef,WvRef,OmvRef] = compute_ref(TththRef,Thth,TthtvRef,Thtv);


Wh_dat(i)   = x(1);
Omh_dat(i)  = x(2);
Thth_dat(i) = x(3);
Wv_dat(i)   = x(4);
Omv_dat(i)  = x(5);
Thtv_dat(i) = x(6);

ref = [WhRef OmhRef TththRef WvRef OmvRef TthtvRef-Thtv0]';

t0 = cputime;
sys = qLPV_TRMS_SS(Wh,Omh,Thth,Wv,Thtv);
mpc = update_mpc_sys_dynamics(mpc,eye(6)+Ts*sys.A,Ts*sys.B,[]);

% Adjust state 6
x_mpc = [Wh;Omh;Thth;Wv;Omv;Thtv-Thtv0];


[u_prev,J,x0] = mpc_solve(x0,x_mpc,u_prev,ref,[],mpc,[],[],[]);
ti(i) = cputime-t0;

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
title('Thth')

figure
plot(Thtv_dat)
grid on
title('Thtv')

figure
plot(uh_dat)
hold on
plot(diff(uh_dat))
grid on
title('uh')

figure
plot(uv_dat)
hold on
plot(diff(uv_dat))
grid on
title('uv')

figure
plot(ti)
title('Compute Time')

