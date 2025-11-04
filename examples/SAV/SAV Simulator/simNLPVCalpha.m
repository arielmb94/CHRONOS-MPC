function x_dot = simNLPVCalpha(w_l,w_r,delta,x,params)
    
m = 1.1937 ;
l_r = 0.1049;
l_f = 0.174-l_r;
t_r = 0.087/2;

C_sigma_2 = params(1);
C_sigma_1 = params(2);
C_sigma_0 = params(3);
C_alpha_f2 = params(4);
C_alpha_f1 = params(5);
C_alpha_f0 = params(6);
C_alpha_r2 = params(7);
C_alpha_r1 = params(8);
C_alpha_r0 = params(9);

Iz = params(13);
R = params(14);
            
vx = x(1);
vy = x(2);
yawrate = x(3);

Fxrr = Fx_funct(w_r,R,vx,C_sigma_2,C_sigma_1,C_sigma_0);
Fxrl = Fx_funct(w_l,R,vx,C_sigma_2,C_sigma_1,C_sigma_0);
Fxr = Fxrr+Fxrl;

Fyf = Fyf_funct(delta,l_f,yawrate,vx,vy,C_alpha_f2,C_alpha_f1,C_alpha_f0);
Fyr = Fyr_funct(l_r,yawrate,vx,vy,C_alpha_r2,C_alpha_r1,C_alpha_r0);

vx_dot = Fxr/m-Fyf*sin(delta)/m+yawrate*vy;
vy_dot = Fyr/m+Fyf*cos(delta)/m-yawrate*vx;
yawrate_dot=(l_f*Fyf*cos(delta)-l_r*Fyr+Fxrr*t_r-Fxrl*t_r)/Iz;

x_dot =  [vx_dot,vy_dot,yawrate_dot]';

end

function Fx = Fx_funct(w,R,vx,C_sigma_2,C_sigma_1,C_sigma_0)
    
    slip_ratio = w*R/vx-1;
    
    C_sigma = C_sigma_2*vx^2+C_sigma_1*vx+C_sigma_0;
    Fx = C_sigma*slip_ratio;

end

function Fyf = Fyf_funct(delta,l_f,yawrate,vx,vy,C_alpha_f2,C_alpha_f1,C_alpha_f0)
        
    slip_angle = delta-atan2(vy+l_f*yawrate,vx);
    
    C_alpha_f = C_alpha_f2*vx^2+C_alpha_f1*vx+C_alpha_f0;
    Fyf = C_alpha_f*slip_angle;
end

function Fyr = Fyr_funct(l_r,yawrate,vx,vy,C_alpha_r2,C_alpha_r1,C_alpha_r0)

    slip_angle = -atan2(vy-l_r*yawrate,vx);
    
    C_alpha_r = C_alpha_r2*vx^2+C_alpha_r1*vx+C_alpha_r0;
    Fyr = C_alpha_r*slip_angle;

end

