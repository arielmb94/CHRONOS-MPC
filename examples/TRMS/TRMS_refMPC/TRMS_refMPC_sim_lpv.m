clear all
TRMS_refMPC_init

%%
tsim = 200;
Sim_samples = tsim/Ts;

TththRef_v = sin(2*pi*(1/57)*(0:Ts:tsim));
TthtvRef_v = -0.6+0.3*sin(2*pi*(1/45)*(0:Ts:tsim));


WhRef = x(1);
WvRef = x(4);

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

[~,OmhRef,~,OmvRef] = compute_ref(TththRef,Thth,TthtvRef,Thtv);

Wh_dat(i)   = x(1);
Omh_dat(i)  = x(2);
Thth_dat(i) = x(3);
Wv_dat(i)   = x(4);
Omv_dat(i)  = x(5);
Thtv_dat(i) = x(6);

ref = [OmhRef TththRef OmvRef TthtvRef-Thtv0]';
x_ref = [WhRef OmhRef TththRef WvRef OmvRef TthtvRef-Thtv0]';

t0 = cputime;
sys = qLPV_TRMS_refMPC_SS(Wh,Omh,Thth,Wv,Thtv);
mpc = update_mpc_sys_dynamics(mpc,eye(6)+Ts*sys.A,Ts*sys.B,[]);

% Adjust state 6
x_mpc = [Wh;Omh;Thth;Wv;Omv;Thtv-Thtv0];


[u_prev,J,x0] = mpc_solve(x0,x_mpc,u_prev,ref,[],mpc,x_ref,[],[]);
ti(i) = cputime-t0;

uh = u_prev(1);
uv = u_prev(2);
WhRef = u_prev(3);
WvRef = u_prev(4);

OmhRef_dat(i) = OmhRef;
OmvRef_dat(i) = OmvRef;

uh_dat(i) = u_prev(1);
uv_dat(i) = u_prev(2);
WhRef_dat(i) = u_prev(3);
WvRef_dat(i) = u_prev(4);


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
plot(ti)
title('Compute Time')

