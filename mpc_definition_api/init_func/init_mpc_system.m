%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = init_mpc_system(mpc,A,B,Bd,C,D,Dd)
%
% Initializes the MPC Discrete-Time dynamical model :
%
%   x+ = A * x + B * u + Bd * d
%
% and the measurement model for the MPC tracking signal:
%
%   y = C * x + D * u + Dd * d
%
% x and u are the state and input vectors, d corresponds to a measured or
% estimated disturbance vector, to be introduced on the appropiate field on
% mpc_solve() during runtime MPC execution.
%
% In:
%   - mpc: CHRONOS mpc structure
%   - A: nx x nx matrix, system matrix
%   - B: nx x nu matrix, input matrix
%   - Bd: nx x nd matrix, disturbance input matrix. If it does not exists,
%   must be set to 0
%   - C: ny x nx matrix, system output matrix
%   - D: ny x nu matrix, output feedtrhough matrix. If it does not exists,
%   must be set to 0
%   - Dd: ny x nd matrix, disturbance output feedtrhough matrix. If it does
%   not exists must be set to 0
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = init_mpc_system(mpc,A,B,Bd,C,D,Dd)

mpc.A = A;
mpc.B = B;
mpc.Bd = Bd;
mpc.C = C;
mpc.D = D;
mpc.Dd = Dd;

mpc.nx = size(A,1);  %number of states
mpc.nu = size(B,2);  %number of control inputs
mpc.nd = size(Bd,2);  %number of disturbance inputs
mpc.ny = size(C,1);  %number of measurements

mpc.Nx = mpc.N*mpc.nx;
mpc.Nu = mpc.N_ctr_hor*mpc.nu;
mpc.Nd = mpc.N*mpc.nd;

if mpc.D == 0
    mpc.Ny = (mpc.N-1)*mpc.ny;
else
    mpc.Ny = mpc.N*mpc.ny;
end

% A equality contraint (b equality constraints depends on x0 and d(k)
mpc.Aeq = genEqualities(A,B,mpc.N,mpc.N_ctr_hor,mpc.Nx,mpc.Nu,mpc.nx,mpc.nu);
mpc.beq = zeros(size(mpc.Aeq,1),1);

end