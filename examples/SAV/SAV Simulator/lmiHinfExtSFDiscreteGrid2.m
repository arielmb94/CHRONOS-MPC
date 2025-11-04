% Author: Ariel Medero
% 
% Description
% Function that computes the Grid Based Hinf SF Discrete-Time LPV controller solving the LMI problem 
% over the BRL.
% To use this function one have to create the generalized plant P st:
%  Xdot    A  Bw  Bu    X
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

function [Kval,Xval,Gval,Yval,CL,gopt] = lmiHinfExtSFDiscreteGrid2(Problem,ncon,percentage,solver)
    
    %%% Define size of the generalized plant
    P = Problem.P{1};
    sizeX = size(P.a,1);
    sizeZ = size(P,1); % number of controlled outputs
    sizeW = size(P,2)-ncon;  % number of input - number of control
    sizeU = ncon;
    
    %%Define the Global Optimization variables: Lyapunov matrix X and Hinf
    %%norm bound gamma
    
    gamma = sdpvar(1,1,'full');
    
    %Create a squared Basis Matrix for each basis element
    %Example: Basis = [1,rho1,rho2] ==> X = X1+rho1*X2+rho*X3
    
    lengthBasis = length(Problem.Basis{1});
    X = CreateGridVariable(sizeX,sizeX,lengthBasis,'symmetric');
    for GridPoint = 1:Problem.TotalPoints
       G{GridPoint} = sdpvar(sizeX,sizeX,'full');
    end
    for GridPoint = 1:Problem.TotalPoints
       Y{GridPoint} = sdpvar(sizeU,sizeX,'full');
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
        Yi = Y{GridPoint};
                
        F = [F,Xi>=0];
        %%% Define LMIs in the Vertex Politope
        for vertex = 1:2^(Problem.NumVaryingParameters)
            
            BasisVariationVertex = Problem.BasisVariation{GridPoint}(vertex,:);
            Xj = GridMatrix(X,BasisVariationVertex);
            
            %%% LMIs definition of the State feedback Hinf problem
            H = Ext_BRL(Xi,Xj,Gi,Yi,gamma,A,Bw,Bu,Cz,Dw,Du,sizeX,sizeX,sizeZ,sizeW);
            
            F = [F,Xj>=0,H>=0];  
        end
    end

    %%% Solve the LMIS
    options = sdpsettings('solver',solver);
    a = optimize(F,gamma,options);

    gopt = double(gamma)
    %%%%%%%
    
    %%% Second iteration to improve numerical aspects with suboptimal gamma
    gnum = gopt*(1+percentage/100);

    %%Define the Global Optimization variables: Lyapunov matrix X and Hinf
    %%norm bound gamma
    
    %Create a squared Basis Matrix for each basis element
    %Example: Basis = [1,rho1,rho2] ==> X = X1+rho1*X2+rho*X3
    clear X Y G
    
    X = CreateGridVariable(sizeX,sizeX,lengthBasis,'symmetric');
    for GridPoint = 1:Problem.TotalPoints
       G{GridPoint} = sdpvar(sizeX,sizeX,'full');
    end
    for GridPoint = 1:Problem.TotalPoints
       Y{GridPoint} = sdpvar(sizeU,sizeX,'full');
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
        Yi = Y{GridPoint};
                
        F = [F,Xi>=0];
        %%% Define LMIs in the Vertex Politope
        for vertex = 1:2^(Problem.NumVaryingParameters)
            
            BasisVariationVertex = Problem.BasisVariation{GridPoint}(vertex,:);
            Xj = GridMatrix(X,BasisVariationVertex);
            
            %%% LMIs definition of the State feedback Hinf problem
            H = Ext_BRL(Xi,Xj,Gi,Yi,gnum,A,Bw,Bu,Cz,Dw,Du,sizeX,sizeX,sizeZ,sizeW);
            
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
       Yval{GridPoint} = double(Y{GridPoint});
    end
    
    %%% Recover the Close-Loop
    for GridPoint = 1:Problem.TotalPoints
        
        P = Problem.P{GridPoint};
        %%% Cut the system
        [A,Bw,Bu,Cz,Dw,Du] = SystemCut(P,sizeX,sizeZ,sizeU,sizeW);
        
        Gi = Gval{GridPoint};
        Yi = Yval{GridPoint};
        
        Ki = Yi/Gi;
        Kval{GridPoint} = Ki; 
        
        CL{GridPoint} = ss(A+Bu*Ki,Bw,Cz+Du*Ki,Dw,Problem.Ts);
    
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

function H = Ext_BRL(X1,X2,G1,Y1,gamma,A,Bw,Bu,Cz,Dw,Du,sizeX1,sizeX2,sizeZ,sizeW)

    M21 = A*G1+Bu*Y1;
    M31 = Cz*G1+Du*Y1;
    M33 = gamma*eye(sizeZ);
    M44 = gamma*eye(sizeW);
    Zero1 = zeros(sizeX1,sizeZ);
    Zero2 = zeros(sizeW,sizeX2);

    H = [G1'+G1-X1 M21' M31' Zero2';
         M21 X2 Zero1 Bw;
         M31 Zero1' M33 Dw;
         Zero2 Bw' Dw' M44];   
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
