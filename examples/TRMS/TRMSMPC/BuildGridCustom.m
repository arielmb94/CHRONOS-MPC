%%%% Step 0: System Parameters definition

%%%% Step 1: Define Grid Points, Parameter Variations and Sampling Time

Problem.NumVaryingParameters = 3;
Wh = [-1.7 -1 -0.5 0.5 1 1.7];
WhRate = 1;
Omh = [-1 1];
OmhRate = 0;
Thth = [-1.5 1.5];
ThthRate = 0;
Wv = [-1.6 -1 -0.5 0.5 1 1.6];
WvRate = 1;
Thtv = [-0.8 -0.45 -0.35 -0.25 0 0.15 0.38];
ThtvRate = 0.8;
Ts = 0.01;                      %s
Problem.Ts = Ts;                

%%%% Step 2: For loop to define the generalazied plant system, basis
%%%% functions and Polytope due to parameter variation at each grid point

%We1 FAN HOR
M = 2;
Ae = 0.01;
f1 = 10; %5
We = ss(tf([1/M 2*pi*f1],[1 2*pi*f1*Ae]));
Wewh = c2d(We,Ts,'Tustin');

M = 2; %HOR ANGLE
Ae = 0.0001;
f1 = 5; %1
We = ss(tf([1/M 2*pi*f1],[1 2*pi*f1*Ae]));
Weh = c2d(We,Ts,'Tustin');

%We2 FAN VERT
M = 2;
Ae = 0.01;
f1 = 8; %3
We = ss(tf([1/M 2*pi*f1],[1 2*pi*f1*Ae]));
Wewv = c2d(We,Ts,'Tustin');

M = 2; %ver ANGLE
Ae = 0.0001;
f2 = 4;  %0.6
We = ss(tf([1/M 2*pi*f2],[1 2*pi*f2*Ae]));
Wev = c2d(We,Ts,'Tustin');

We = [Wewh 0 0 0;0 Weh 0 0;0 0 Wewh 0; 0 0 0 Wev];
We.u = 'e';
We.y = 'ze';

%Wu1  HOR
f3 = 20; %10
Mu = 2;
Wu = ss(tf([1 2*pi*f3/Mu], [0.001 2*pi*f3]));
Wuh = c2d(Wu,Ts,'Tustin');

%Wu2   VERT
f4 = 15; %6
Mu = 1.5;
Wu = ss(tf([1 2*pi*f4/Mu], [0.001 2*pi*f4]));
Wuv = c2d(Wu,Ts,'Tustin');

Wu = [Wuh 0;0 Wuv];
Wu.u = 'u';
Wu.y = 'zu';

%%
tic
GridPoint = 0;
for WhPoint = 1:length(Wh)
    
%wh value at the grid point
wh = Wh(WhPoint);
    
for OmhPoint = 1:length(Omh)
        
%wh value at the grid point
omh = Omh(OmhPoint);  

for ThthPoint = 1:length(Thth)

%thth value at the grid point
thth = Thth(ThthPoint);

for WvPoint = 1:length(Wv)
  
%wv value at the grid point    
wv = Wv(WvPoint);   

for ThtvPoint = 1:length(Thtv)
    
%thtv value at the grid point
thtv = Thtv(ThtvPoint);

%System at GridPoint
G = qLPV_TRMS_SS(wh,omh,thth,wv,thtv);
Gd = c2d(G,Ts,'zoh');
Gd.u = 'u';
Gd.y = 'y';

%Increase Grid Point Counter
GridPoint = GridPoint+1;

%Generalized Plant at Grid Point
S = sumblk('e = r - y',4);
P = connect(Gd,We,Wu,S,{'r';'u'},{'ze';'zu'});
Problem.P{GridPoint} = P; % Generalized Plant at GridPoint

%Basis Function
Problem.Basis{GridPoint} = [1,wh,wv,thtv,cos(thtv),sin(thtv)];

%Parameter Variation Politope Basis
whUp = wh+Ts*WhRate;
whDown = wh-Ts*WhRate;
wvUp = wv+Ts*WvRate;
wvDown = wv-Ts*WvRate;
thtvUp = thtv+Ts*ThtvRate;
thtvDown = thtv-Ts*ThtvRate;

Problem.BasisVariation{GridPoint} = [1,whUp,wvUp,thtvUp,cos(thtvUp),sin(thtvUp);
                                     1,whUp,wvUp,thtvDown,cos(thtvDown),sin(thtvDown);
                                     1,whUp,wvDown,thtvUp,cos(thtvUp),sin(thtvUp);
                                     1,whUp,wvDown,thtvDown,cos(thtvDown),sin(thtvDown);
                                     1,whDown,wvUp,thtvUp,cos(thtvUp),sin(thtvUp);
                                     1,whDown,wvUp,thtvDown,cos(thtvDown),sin(thtvDown);
                                     1,whDown,wvDown,thtvUp,cos(thtvUp),sin(thtvUp);
                                     1,whDown,wvDown,thtvDown,cos(thtvDown),sin(thtvDown)];
                                     
end
    
end
    
end
        
end
        
end

Problem.TotalPoints = GridPoint;
toc

 %%
tic
[K,X,gamma] = lmiHinfSFDiscreteGrid_PL(Problem,2,15,'mosek')
toc

Z.Kval = K;
%% Extended BRL with G
tic
[K,Xval,Gval,CL,gopt] = lmiHinfExtSFDiscreteGrid_PL_V2(Problem,2,15,'mosek')
toc
Z.Kval = K;
 %% FL Approach
tic
[K,X,sigma,gamma,CL] = lmiHinfSFDiscreteGrid_FL(Problem,2,15,'mosek',1)
toc

Z.Kval = K;
%%
%X = Xcont;

Keval = K{1}+K{2}*wh+K{3}*wv+K{4}*thtv+K{5}*cos(thtv)+K{6}*sin(thtv)
Xeval = X{1}+X{2}*wh+X{3}*wv+X{4}*thtv+X{5}*cos(thtv)+X{6}*sin(thtv);
eig(Xeval)

% Keval = K{1}+wh*K{2}+K{3}*(sign(wh)*wh^2)+K{4}*omh+K{5}*thth+K{6}*wv+K{7}*(sign(wv)*wv^2)+K{8}*thtv+K{9}*cos(thtv)+K{10}*sin(thtv)
% Xeval = X{1}+wh*X{2}+X{3}*(sign(wh)*wh^2)+X{4}*omh+X{5}*thth+X{6}*wv+X{7}*(sign(wv)*wv^2)+X{8}*thtv+X{9}*cos(thtv)+X{10}*sin(thtv);
% eig(Xeval)

