nx = 6; %[il wl ir wr]
nw = 5;
nu = 2;
nz = 6;
ny = 3;

A = [-Ra_l/La_l -K_l/La_l 0 0 0 0;
     K_l/J_l -F_l/J_l 0 0 0 0;
     0 0 -Ra_r/La_r -K_r/La_r 0 0;
     0 0 K_r/J_r -F_r/J_r 0 0;
     zeros(2,6)];    
B = [1/La_l 0;0 0;0 1/La_r;0 0;zeros(2,2)];
Bw = [zeros(1,nw);-0.01/J_l zeros(1,nw-1);zeros(1,nw);0 -0.01/J_r zeros(1,nw-2);zeros(2,nw)];
C = [0 1 0 0 1 0;
     0 0 0 1 0 1;
     1 0 1 0 0 0];
Cobs = C;
Dw = [zeros(ny,nw-3) diag([5000 5000 10])];%10 10 5
D = zeros(3,2);
sys = ss(A,[Bw B],C,[Dw D]);
ts = 0.02;
sysd = c2d(sys,ts,'Tustin');




%%

options = sdpsettings('solver','sdpt3');

alpha = 0.95;%1;

A = sysd.A;
Aobs = A;
Bw = sysd.B(:,1:nw);
Bu = sysd.B(:,nw+1:nw+nu);
Buobs = Bu;
Cz = diag([1 1 1 1 100 100]);
Dzw = zeros(nz,nw);

X = sdpvar(nx,nx);
G = sdpvar(nx,nx,'full');
Y = sdpvar(nx,ny,'full');
gamma = sdpvar;

H = [G+G'-X G'*A-Y*C G'*Bw-Y*Dw zeros(nx,nz);
     (G'*A-Y*C)' alpha*X zeros(nx,nw) Cz';
     (G'*Bw-Y*Dw)' zeros(nw,nx) gamma*eye(nw) Dzw';
     zeros(nz,nx) Cz Dzw gamma*eye(nz)];
  
Fopt = [X>=0,H>=0,G+G'-X>=0];
%optimize(Fopt,gamma,options);
gopt = double(gamma);
gnum = gopt*1.15;

X = sdpvar(nx,nx);
G = sdpvar(nx,nx,'full');
Y = sdpvar(nx,ny,'full');

H = [G+G'-X G'*A-Y*C G'*Bw-Y*Dw zeros(nx,nz);
     (G'*A-Y*C)' alpha*X zeros(nx,nw) Cz';
     (G'*Bw-Y*Dw)' zeros(nw,nx) gnum*eye(nw) Dzw';
     zeros(nz,nx) Cz Dzw gnum*eye(nz)];
  
Fopt = [X>=0,H>=0];

%optimize(Fopt,[],options);

Xv = double(X);
Gv = double(G);
Yv = double(Y);

Lobs = inv(Gv)'*Yv
gopt

%%
A = [-Ra_l/La_l -K_l/La_l 0 0;
     K_l/J_l -F_l/J_l 0 0;
     0 0 -Ra_r/La_r -K_r/La_r;
     0 0 K_r/J_r -F_r/J_r];    
B = [1/La_l 0;0 0;0 1/La_r;0 0];
Bw = [0 0;-1/J_l 0;0 0;0 -1/J_r];
C = [0 1 0 0;
     0 0 0 1;
     1 0 1 0];
Gmotor = ss(A,[Bw B],C,0);
Gmotor.InputDelay = [0 0 tau_m tau_m];
