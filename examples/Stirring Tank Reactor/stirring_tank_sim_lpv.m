% Load init script defining the mpc problem
stirring_tank_init

clear c_dat v_dat uk rk
ck = x_prev(1);
vk = x_prev(2);


t1 = cputime;
for i = 1:2400

if i < 900
    ref_c = 0.27 + (0.65-0.27)*i/900;
elseif i >= 900 && i<1800
    ref_c = 0.65;
else
    ref_c = 0.65 - (0.65-0.3)*(i-1800)/600;
end

ref_v = - M/(log(1/(theta_f*k*ref_c)*(1-ref_c)));

% Update LPV model
A_lpv = eye(2)+Ts*[-1/theta_f-k*exp(-M/vk) -k*ck*M*exp(-M/vk)/(vk^2);
     k*exp(-M/vk) -1/theta_f];

B_lpv = Ts*[0; -alpha*(vk-xc)];

Bd_lpv = Ts*[1/theta_f k*ck*M*exp(-M/vk)/(vk^2); xf/theta_f 0];

mpc = update_mpc_sys_dynamics(mpc,A_lpv,B_lpv,Bd_lpv);

d = [1;vk];

xref = [ref_c;ref_v];

[u_prev,J,x0] = mpc_solve(x0,x_prev,u_prev,xref,d,mpc,[],[],[]);

ck = ck + Ts*((1-ck)/theta_f - k*ck*exp(-M/vk));
vk = vk + Ts*((xf-vk)/theta_f + k*ck*exp(-M/vk)-alpha*u_prev*(vk-xc));

x_prev = [ck;vk];

c_dat(:,i) = ck;
v_dat(:,i) = vk;
Jk(i) = J;
uk(i) = u_prev;
rk(i) = ref_c;

end
t2 = cputime-t1;
t2/2400
2400/t2

close all

figure
plot(Ts*(1:i),rk,'r',Ts*(1:i),c_dat,'g',Ts*(1:i),v_dat,'b')
grid on
ylim([0 1])
legend('c ref','ck','vk')

figure
plot(Ts*(1:i),uk,Ts*(1:i-1),diff(uk))
grid on
legend('U','Delta U')