clear all
TRMS_two_lvl_init

%%
tsim = 200;
Sim_samples = tsim/Ts;

TththRef_v = sin(2*pi*(1/57)*(0:Ts:tsim));
TthtvRef_v = -0.6+0.3*sin(2*pi*(1/45)*(0:Ts:tsim));

% TththRef_v = 0.3*ones(1,Sim_samples);
% TthtvRef_v = -0.4*ones(1,Sim_samples);


%%
clear Wh_dat Omh_dat Thth_dat Wv_dat Omv_dat Thtv_dat uh_dat uv_dat ti ...
    WhRef_dat WvRef_dat OmhRef_dat OmvRef_dat

for i = 1:Sim_samples

Wh   = x(1);
Omh  = x(2);
Thth = x(3);
Wv   = x(4);
Omv  = x(5);
Thtv = x(6);

TthtvRef = TthtvRef_v(i);
TththRef = TththRef_v(i);

OmhRef = (TththRef-Thth)/3;
OmvRef = (TthtvRef-Thtv)/5;

Wh_dat(i)   = x(1);
Omh_dat(i)  = x(2);
Thth_dat(i) = x(3);
Wv_dat(i)   = x(4);
Omv_dat(i)  = x(5);
Thtv_dat(i) = x(6);

t0 = cputime;

[A,B,Bd,Ah,Bh,Av,Bv] = qLPV_TRMS_two_lvl_SS(Wh,Omh,Thth,Wv,Thtv);

mpc = update_mpc_sys_dynamics(mpc,eye(4)+Ts*A,Ts*B,Ts*Bd);

mpc_h = update_mpc_sys_dynamics(mpc_h,1+Ts*Ah,Ts*Bh,[]);

mpc_v = update_mpc_sys_dynamics(mpc_v,1+Ts*Av,Ts*Bv,[]);

% Adjust state 6
x_mpc = [Omh;Thth;Omv;Thtv-Thtv0];
ref = [OmhRef TththRef OmvRef TthtvRef-Thtv0]';

[u_prev,J,x0] = mpc_solve(x0,x_mpc,u_prev,ref,uv,mpc,[],[],[]);

WhRef = u_prev(1);
WvRef = u_prev(2);

[uh,J,x0_h] = mpc_solve(x0_h,Wh,uh,WhRef,[],mpc_h,[],[],[]);
[uv,J,x0_v] = mpc_solve(x0_v,Wv,uv,WvRef,[],mpc_v,[],[],[]);
ti(i) = cputime-t0;

% WhRef = u_prev(1);
% WvRef = u_prev(2);

% OmhRef_dat(i) = OmhRef;
% OmvRef_dat(i) = OmvRef;
% 
WhRef_dat(i) = WhRef;
WvRef_dat(i) = WvRef;

uh_dat(i) = uh;
uv_dat(i) = uv;

dt_x = TRMS(Wh,Omh,Thth,Wv,Omv,Thtv,uh,uv);
x = x + Ts*dt_x;

end

1/mean(ti)
%%
close all
plot(TththRef_v)
hold on
plot(Thth_dat)
grid on
title('Thth')

figure
plot(TthtvRef_v)
hold on
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
plot(WhRef_dat)
hold on
plot(diff(WhRef_dat))
grid on
title('WhRef')

figure
plot(WvRef_dat)
hold on
plot(diff(WvRef_dat))
grid on
title('WvRef')

figure
plot(ti)
title('Compute Time')

