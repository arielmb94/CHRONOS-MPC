% Author: Ariel Medero
% 
% Description:
% Function that computes the Hinf SF Discrete time LTI controller solving the LMI problem 
% over the BRL.
% To use this function one have to create the generalized plant P st:
%  Xk+1    A  Bw  Bu    X
%   Z   =  Cz Dzw Dzu   W
%   Y      Cy Dyw  0    U
% 
% Input
%  P  : LTI generalized plant
%  ncon   : number of control signals
%  percentage : percentage added to the gamma optimal to conditionate the controller
%  solver : name of the solver used ('dsdp','sedumi'...)
% 
% Output
%  K  : Controller (state space)
%  gopt   : optimal gamma
%
% Use:
% [K,gopt] = lmiHinfSFDiscrete(P,ncon,percentage,solver)

function [K,gopt,gval,Pcl,CL,Xval,Yval] = lmiHinfDOFDiscrete_PL(P,ncon,nmeas,sp,percentage,solver)

    %%% Size of the generalized plant
    sizeX = size(P.a,1);
    sizeZ = size(P,1)-nmeas; % number of controlled output - number of measure
    sizeY = nmeas;
    sizeW = size(P,2)-ncon;  % number of input - number of control
    sizeU = ncon;

    %%% Cut the system
    [A,Bw,Bu,Cz,Dzw,Dzu,Cy,Dyw] = SystemCut(P,sizeX,sizeZ,sizeY,sizeU,sizeW);

    %%% Create variables matrix
    X = sdpvar(sizeX,sizeX); 
    Y = sdpvar(sizeX,sizeX);
    gamma = sdpvar(1,1,'full');

    %%% LMIs definition of the State feedback Hinf problem
    [H1,H2,H3] = BRL_PL(X,X,Y,Y,gamma,A,Bw,Bu,Cz,Dzw,Dzu,Cy,Dyw,sizeX,sizeZ,sizeW);

    F = [H1<=0,H2<=0,H3>=0];
    %%% Solve the LMIS
    options = sdpsettings('solver',solver);
    a = optimize(F,gamma,options);

    gopt = double(gamma);
    
    %%% Second iteration to improve numerical aspects with suboptimal gamma
    gnum = gopt*(1+percentage/100);
    
    %%% Create variables matrix
    clear X Y
    X = sdpvar(sizeX,sizeX); 
    Y = sdpvar(sizeX,sizeX);

    %%% LMIs definition of the Hinf problem
    [H1,H2,H3] = BRL_PL(X,X,Y,Y,gnum,A,Bw,Bu,Cz,Dzw,Dzu,Cy,Dyw,sizeX,sizeZ,sizeW);

    F = [H1<=0,H2<=0,H3>=0];
    %%% Solve the LMIS
    options = sdpsettings('solver',solver);
    a = optimize(F,[],options);

    %%% Optimization over the BRL in CL to obtain the controller
    Xval = double(X);
    Yval = double(Y);
    
    %%% Back to the solution
    [u,s,v] = svd(eye(sizeX) - Xval*Yval);
    M = u*s.^0.5;
    N = v*s.^0.5;
    nro = size(s,1);
    Pcl = [Yval eye(sizeX);N' zeros(nro,sizeX)]*inv([eye(sizeX) Xval;zeros(nro,sizeX) M']);
    Pcl = triu(Pcl)+triu(Pcl,1)';
%     Xclinv = [R eye(sizeX); M' zeros(nro,sizeX)] / [eye(sizeX) S; zeros(nro,sizeX) N']
%     [U,S,V] = svd(eye(sizeX) - Xval*Yval);
%     R       = chol(S); % R'*R=S
%     M       = U*R';
%     N       = (R*V')';
    
%     Pcl = [Yval eye(sizeX);N' zeros(sizeX)]/[eye(sizeX) Xval;zeros(sizeX) M'];
%     max_num_error = max(max(Pcl-Pcl'));
%     power=ceil(log10(max_num_error)-1)+2;
%     
%     Pcl = round(Pcl,abs(power));
    
    Ak = sdpvar(sizeX,sizeX,'full');
    Bk = sdpvar(sizeX,sizeY,'full');
    Ck = sdpvar(sizeU,sizeX,'full');
    if sp
        Dk = zeros(sizeU,sizeY);
    else
        Dk = sdpvar(sizeU,sizeY,'full');
    end
    
    gamma = sdpvar(1,1,'full');
    H = BRL_CL(Pcl,Pcl,gamma,Ak,Bk,Ck,Dk,A,Bw,Bu,Cz,Dzw,Dzu,Cy,Dyw,sizeX,sizeZ,sizeW);
    F = H<=0;

    options = sdpsettings('solver',solver);
    a = optimize(F,gamma,options)

    gval = value(gamma);
    
    Akval = double(Ak);
    Bkval = double(Bk);
    Ckval = double(Ck);
    Dkval = double(Dk);
    K = ss(Akval,Bkval,Ckval,Dkval,P.Ts);

    CL = lft(P,K);

end

function [A,Bw,Bu,Cz,Dzw,Dzu,Cy,Dyw] = SystemCut(P,sizeX,sizeZ,sizeY,sizeU,sizeW)
    A   = P.a(1:sizeX,1:sizeX);
    Bw  = P.b(1:sizeX,1:sizeW);
    Bu  = P.b(1:sizeX,sizeW+1:sizeW+sizeU);
    Cz  = P.c(1:sizeZ,1:sizeX);
    Dzw = P.d(1:sizeZ,1:sizeW);
    Dzu = P.d(1:sizeZ,sizeW+1:sizeW+sizeU);
    Cy  = P.c(sizeZ+1:sizeZ+sizeY,1:sizeX);
    Dyw = P.d(sizeZ+1:sizeZ+sizeY,1:sizeW);
    D22 = P.d(sizeZ+1:sizeZ+sizeY,sizeW+1:sizeW+sizeU);
end

function [H1,H2,H3] = BRL_PL(X1,X2,Y1,Y2,gamma,A,Bw,Bu,Cz,Dzw,Dzu,Cy,Dyw,sizeX,sizeZ,sizeW)

    %First Projector
    rW = rank([Bu' Dzu']);
    Nw = null([Bu' Dzu']);
    ProjW = [Nw zeros(sizeX+sizeZ,sizeW);
             zeros(sizeW,sizeX+sizeZ-rW) eye(sizeW)];
    
    H = [A*X2*A'-X1 A*X2*Cz' Bw;
         Cz*X2*A' Cz*X2*Cz'-gamma*eye(sizeZ) Dzw;
         Bw' Dzw' -gamma*eye(sizeW)];
     
    H1 = ProjW'*H*ProjW;
         
    %Second Projector
    rY = rank([Cy Dyw]);
    Ny = null([Cy Dyw]);
    ProjY = [Ny zeros(sizeX+sizeW,sizeZ);
             zeros(sizeZ,sizeX+sizeW-rY) eye(sizeZ)];
    
    H = [A'*Y2*A-Y1 A'*Y2*Bw Cz';
         Bw'*Y2*A Bw'*Y2*Bw-gamma*eye(sizeW) Dzw';
         Cz Dzw -gamma*eye(sizeZ)];   

    H2 = ProjY'*H*ProjY;
    
    H3 = [X2 eye(sizeX);eye(sizeX) Y2];
    
end

function H = BRL_CL(Pcl1,Pcl2,gamma,Ak,Bk,Ck,Dk,A,Bw,Bu,Cz,Dzw,Dzu,Cy,Dyw,sizeX,sizeZ,sizeW)

    Acl = [A+Bu*Dk*Cy Bu*Ck;Bk*Cy Ak];
    Bcl = [Bw+Bu*Dk*Dyw;Bk*Dyw];
    Ccl = [Cz+Dzu*Dk*Cy Dzu*Ck];
    Dcl = Dzw+Dzu*Dk*Dyw;
     
    H =[-inv(Pcl2),Acl,Bcl,zeros(sizeX*2,sizeZ);
        Acl',-Pcl1,zeros(sizeX*2,sizeW),Ccl';
        Bcl',zeros(sizeW,sizeX*2),-gamma*eye(sizeW),Dcl';
        zeros(sizeZ,sizeX*2),Ccl,Dcl,-gamma*eye(sizeZ)]
    
end