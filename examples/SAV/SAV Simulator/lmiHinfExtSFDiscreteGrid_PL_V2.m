% Author: Ariel Medero
% 
% Description
% Function that computes the Grid Based Hinf SF Discrete-Time LPV controller solving the LMI problem 
% over the BRL.
% To use this function one have to create the generalized plant P st:
%  Xk+1    A  Bw  Bu    X
%   Z   =  Cz Dw  Du    W
% 
% Input
%  P  : LTI generalized plant
%  ncon   : number of control signals
%  percentage : percentage added to the gamma optimal to conditionate the controller
%  solver : name of the solver used ('dsdp','sedumi'...)
% 
% Output
%  K  : Controller (state space)
%  CL : closed loop (state space)
%  gopt   : optimal gamma
%
% [K,CL,gopt] = lmiHinfSF(P,nstate,ncon,percentage,solver)

function [Kval,Xval,Gval,CL,gopt] = lmiHinfExtSFDiscreteGrid_PL_V2(Problem,ncon,percentage,solver)
    
    %%% Define size of the generalized plant
    P = Problem.P{1}
    sizeX = size(P.a,1)
    sizeZ = size(P,1) % number of controlled outputs
    sizeW = size(P,2)-ncon  % number of input - number of control
    sizeU = ncon
    
    %%Define the Global Optimization variables: Lyapunov matrix X and Hinf
    %%norm bound gamma
    
    gamma2 = sdpvar(1,1,'full');
    
    %Create a squared Basis Matrix for each basis element
    %Example: Basis = [1,rho1,rho2] ==> X = X1+rho1*X2+rho*X3
    
    lengthBasis = length(Problem.Basis{1});
    X = CreateGridVariable(sizeX,sizeX,lengthBasis,'symmetric');
    
    %G unique for each GridPoint
    for GridPoint = 1:Problem.TotalPoints
        G{GridPoint} = sdpvar(sizeX,sizeX,'full');        
    end
    
    %Create a for loop to define the LMIs at each Gridpoint
    F = [];
    for GridPoint = 1:Problem.TotalPoints

        P = Problem.P{GridPoint};
        %%% Cut the system
        [A,Bw,Bu,Cz,Dw,Du] = SystemCut(P,sizeX,sizeZ,sizeU,sizeW);
        
        %%% Lyapunov Matrix at gridpoints
        Xi = GridMatrix(X,Problem.Basis{GridPoint});
        Gi = G{GridPoint};
        
        F = [F,Xi>=0];
        %%% Define LMIs in the Vertex Politope
        for vertex = 1:2^(Problem.NumVaryingParameters)
            
            BasisVariationVertex = Problem.BasisVariation{GridPoint}(vertex,:);
            Xj = GridMatrix(X,BasisVariationVertex);
            
            %%% LMIs definition of the State feedback Hinf problem
            H = BRL_PL(Xi,Xj,Gi,gamma2,A,Bw,Bu,Cz,Dw,Du,sizeX,sizeZ,sizeW);
            
            F = [F,Xj>=0,H>=0];           
        end
    end

    %%% Solve the LMIS
    options = sdpsettings('solver',solver);
   % options.sdpt3.maxit=200;
    a = optimize(F,gamma2,options);

    gopt = sqrt(value(gamma2));
    
    %%%%%%%
    
    %%% Second iteration to improve numerical aspects with suboptimal gamma
    gnum = gopt*(1+percentage/100);
    gnum2 = gnum^2;
    %%Define the Global Optimization variables: Lyapunov matrix X and Hinf
    %%norm bound gamma
    
    %Create a squared Basis Matrix for each basis element
    %Example: Basis = [1,rho1,rho2] ==> X = X1+rho1*X2+rho*X3
    clear X G
    X = CreateGridVariable(sizeX,sizeX,lengthBasis,'symmetric');
   
        %G unique for each GridPoint
    for GridPoint = 1:Problem.TotalPoints
        G{GridPoint} = sdpvar(sizeX,sizeX,'full');        
    end
   
    %Create a for loop to define the LMIs at each Gridpoint
    F = [];
    for GridPoint = 1:Problem.TotalPoints

        P = Problem.P{GridPoint};
        %%% Cut the system
        [A,Bw,Bu,Cz,Dw,Du] = SystemCut(P,sizeX,sizeZ,sizeU,sizeW);
        
        %%% Lyapunov Matrix at gridpoints
        Xi = GridMatrix(X,Problem.Basis{GridPoint});
        Gi = G{GridPoint};

        F = [F,Xi>=0];
        %%% Define LMIs in the Vertex Politope
        for vertex = 1:2^(Problem.NumVaryingParameters)
            
            BasisVariationVertex = Problem.BasisVariation{GridPoint}(vertex,:);
            Xj = GridMatrix(X,BasisVariationVertex);
            
            %%% LMIs definition of the State feedback Hinf problem
            H = BRL_PL(Xi,Xj,Gi,gnum2,A,Bw,Bu,Cz,Dw,Du,sizeX,sizeZ,sizeW);
            
            F = [F,Xj>=0,H>=0];  
        end
    end

    %%% Solve the LMIS
    options = sdpsettings('solver',solver);
    a = optimize(F,[],options);

    %%% Recover the Lyapunov Matrix solution
    for i = 1:length(Problem.Basis{1})
        Xval{i} = double(X{i});
    end
    
    for GridPoint = 1:Problem.TotalPoints
        Gval{GridPoint} = double(G{GridPoint});        
    end
    
    %%%%%%%
    
    %%% Optimization over the BRL in CL to obtain the controller
    K = CreateGridVariable(sizeU,sizeX,lengthBasis,'full');
    gamma2 = sdpvar(1,1,'full');

    F = [];
    for GridPoint = 1:Problem.TotalPoints

        P = Problem.P{GridPoint};
        %%% Cut the system
        [A,Bw,Bu,Cz,Dw,Du] = SystemCut(P,sizeX,sizeZ,sizeU,sizeW);
        
        %%% Lyapunov Matrix at gridpoint
        Xi = GridMatrix(Xval,Problem.Basis{GridPoint});
        Gi = Gval{GridPoint};
        Ki = GridMatrix(K,Problem.Basis{GridPoint});
        
        %%% Define LMIs in the Vertex Politope
        for vertex = 1:2^(Problem.NumVaryingParameters)
            
            BasisVariationVertex = Problem.BasisVariation{GridPoint}(vertex,:);
            Xj = GridMatrix(Xval,BasisVariationVertex);
            
            %%% LMIs definition of the State feedback Hinf problem
            H = BRL_CL(Xi,Xj,Gi,Ki,gamma2,A,Bw,Bu,Cz,Dw,Du,sizeX,sizeZ,sizeW);
            
            F = [F,H>=0];  
        end
    end

    %%% Solve the LMIS
    options = sdpsettings('solver',solver);
    a = optimize(F,gamma2,options);
    
    %%% Recover the Controller solution
    for i = 1:length(Problem.Basis{1})
        Kval{i} = double(K{i});
    end
    
    %%% Recover the Close-Loop at each grid point
    for GridPoint = 1:Problem.TotalPoints
        
       P = Problem.P{GridPoint};
       %%% Cut the system
       [A,Bw,Bu,Cz,Dw,Du] = SystemCut(P,sizeX,sizeZ,sizeU,sizeW);
       
       Ki = GridMatrix(Kval,Problem.Basis{GridPoint});
       
       CL{GridPoint} = ss(A+Bu*Ki,Bw,Cz+Du*Ki,Dw,P.Ts);
        
    end
    
end

function [A,Bw,Bu,Cz,Dw,Du] = SystemCut(P,sizeX,sizeZ,sizeU,sizeW)
    A   = P.a(1:sizeX,1:sizeX);
    Bw  = P.b(1:sizeX,1:sizeW);
    Bu  = P.b(1:sizeX,sizeW+1:sizeW+sizeU);
    Cz  = P.c(1:sizeZ,1:sizeX);
    Dw = P.d(1:sizeZ,1:sizeW);
    Du = P.d(1:sizeZ,sizeW+1:sizeW+sizeU);
end

function H1 = BRL_PL(X1,X2,G,gamma2,A,Bw,Bu,Cz,Dw,Du,sizeX,sizeZ,sizeW)

    %First Projector
    rQ = rank([Bu' Du']);
    NQ = null([Bu' Du']);
    ProjQ = [eye(sizeX) zeros(sizeX,sizeX+sizeZ-rQ) zeros(sizeX,sizeW);
             zeros(sizeX+sizeZ,sizeX) NQ zeros(sizeX+sizeZ,sizeW);
             zeros(sizeW,sizeX) zeros(sizeW,sizeX+sizeZ-rQ) eye(sizeW)];
             
    H = [G+G'-X2 G'*A' G'*Cz' zeros(sizeX,sizeW);
         (G'*A')' X1 zeros(sizeX,sizeZ) Bw;
         (G'*Cz')' zeros(sizeZ,sizeX) gamma2*eye(sizeZ) Dw;
         zeros(sizeW,sizeX) Bw' Dw' eye(sizeW)];
     
    H1 = ProjQ'*H*ProjQ;
    
end

function H = BRL_CL(X1,X2,G,K,gamma2,A,Bw,Bu,Cz,Dw,Du,sizeX,sizeZ,sizeW)

    H =[G+G'-X2 G'*A'+G'*K'*Bu' G'*Cz'+G'*K'*Du' zeros(sizeX,sizeW);
        (G'*A'+G'*K'*Bu')' X1 zeros(sizeX,sizeZ) Bw;
        (G'*Cz'+G'*K'*Du')' zeros(sizeZ,sizeX) gamma2*eye(sizeZ) Dw;
        zeros(sizeW,sizeX) Bw' Dw' eye(sizeW)];   
    
end

function Mi = GridMatrix(M,BasisVec)
%%% This Function evaluate a Parameter Dependent Matrix at a GridPoint
    [dimA,dimB] = size(M{1});
    Mi = zeros(dimA,dimB);
    for i = 1:length(BasisVec)
        Mi = Mi+BasisVec(i)*M{i}; 
    end
end

function M = CreateGridVariable(dimA,dimB,lengthBasis,type)
%%% This Function creates a Parameter Dependent Matrix
%%% Type is a string indicating which type of Matrix is: 'full' or
%%% 'symmetric'

    for i = 1:lengthBasis
        M{i} = sdpvar(dimA,dimB,type); 
    end
end
