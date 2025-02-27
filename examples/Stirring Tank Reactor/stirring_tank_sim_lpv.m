% Load init script defining the mpc problem
stirring_tank_init

clear c_dat v_dat uk rk
ck = x_prev(1);
vk = x_prev(2);

r_c = 0.65;
r_v = 0.53;

f_c = ck;
f_v = vk;
tau = 15;

t1 = cputime;
for i = 1:2400

    f_c = f_c + Ts*(-f_c/tau+r_c/tau);
    f_v = f_v + Ts*(-f_v/tau+r_v/tau);
    
    % Update LPV model
    A_lpv = eye(2)+Ts*[-1/theta_f-k*exp(-M/vk) -k*ck*M*exp(-M/vk)/(vk^2);
         k*exp(-M/vk) -1/theta_f];
    
    B_lpv = Ts*[0; -alpha*(vk-xc)];
    
    Bd_lpv = Ts*[1/theta_f k*ck*M*exp(-M/vk)/(vk^2); xf/theta_f 0];
    
    mpc = update_mpc_sys_dynamics(mpc,A_lpv,B_lpv,Bd_lpv);
    
    d = [1;vk];
    
    xref = [f_c;f_v];
    
    [u_prev,J,x0,fallback_control] = mpc_solve(x0,x_prev,u_prev,xref,d,mpc,[],[],[]);
    if fallback_control
        if i > 1
            if Jk(i - 1) < J
                u_prev = uk(i - 1);            
            end
            disp("Fallback control is used!")
        end
    end
    
    ck = ck + Ts*((1-ck)/theta_f - k*ck*exp(-M/vk));
    vk = vk + Ts*((xf-vk)/theta_f + k*ck*exp(-M/vk)-alpha*u_prev*(vk-xc));
    
    x_prev = [ck;vk];
    
    c_dat(:,i) = ck;
    v_dat(:,i) = vk;
    Jk(i) = J;
    uk(i) = u_prev;
    rk(i) = f_c;

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