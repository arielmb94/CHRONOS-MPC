function [ TheMPC] = cMPC(A,B,C,Np,nx,ny,nu,Qy,R,xmin,xmax,umin,umax)

%SDPVARS
x          = sdpvar(nx,Np);
y          = sdpvar(ny,Np-1); 
u          = sdpvar(nu,Np-1); 


%PARAMETERS
yref        = sdpvar(ny,1);
uref        = sdpvar(nu,1); 


% Objective function and constraints   
    J = 0;   
    st = []; 
    
for k = 1:Np-1   
    J = J +(y(:,k)-yref)'*Qy*(y(:,k)-yref)+(u(:,k)-uref)'*R*(u(:,k)-uref);    
    st = [st, x(:,k+1) == A*x(:,k) + B*u(:,k)];
    st = [st, y(:,k)   == C*x(:,k)];
    st = [st, xmin <= x(:,k) <= xmax];
    st = [st, umin <= u(:,k) <= umax];
end
    st = [st, xmin <= x(:,Np) <= xmax];

options = sdpsettings('verbose',1,'solver','SEDUMI');

%options = sdpsettings('verbose',1,'solver','fmincon','fmincon.Algorithm','interior-point','fmincon.SubproblemAlgorithm','cg','fmincon.Display','off');

Parameters = {x(:,1) yref uref};

solutions_out = {u(:,1)};
TheMPC = optimizer(st,J,options,Parameters,solutions_out); 

end

