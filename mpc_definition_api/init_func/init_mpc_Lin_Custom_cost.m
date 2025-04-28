%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Adds quadratic and linear penalties on a custom user defined signal z:
%
%   J += z'*Qz*z + qz*z
%
% The signal z is defined as:
%
%   z = Cz * x + Dz * u + Ddz * dz
%
% The user must define the signal z by selecting appropiate values for the
% matrices Cz, Dz, Ddz.
%
% In:
%   - mpc: CHRONOS mpc structure
%   - Cz: nz x nx matrix, states output matrix
%   - Dz: nz x nu matrix, input feedtrhough matrix
%   - Ddz: nz x ndi matrix, disturbance feedtrhough matrix
%   - Qz (optional): nz x nz square matrix, weights for the quadratic
%   penalty term on the user defined signal z
%   - qz (optional): nz column vector, weights for the linear penalty term
%   on the user defined signal z. Use the linear penalty term only in the
%   case that z takes only positive values
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
% Example:
% We want to penalize the variation of the MPC control action with respect 
% a given value u_star, e.g. we want to minimize z = u - u_star.
% To achieve this we define z by selecting:
%   - Cz = [0 0 ... 0]
%   - Dz = [1]
%   - Ddz = [-1]
% which corresponds to z = [0 0 ... 0] * x + [1] * u + [-1] * di
% The "disturbance" term on z corresponds to u_star, to be introduced on 
% the appropiate field on mpc_solve() during runtime MPC execution.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = init_mpc_Lin_Custom_cost(mpc,Cz,Dz,Ddz,Qz,qz)
arguments
    mpc
    Cz
    Dz
    Ddz
    Qz = [];
    qz = [];
end

mpc.Qz = Qz;
mpc.qz = qz;

mpc.Cz = Cz;
mpc.Dz = Dz;
mpc.Ddz = Ddz;

mpc.ndz = size(Ddz,2);  %number of disturbance inputs to performance cost
mpc.nz = size(Cz,1);  %number of performances

if mpc.Dz == 0
    mpc.Nz = (mpc.N-1)*mpc.nz;
else
    mpc.Nz = mpc.N*mpc.nz;
end

mpc.Ndz = mpc.N*mpc.ndz;

% Quadratic Cost Term Gradient and Hessian Computation
if ~isempty(Qz)
    
    if isempty(mpc.hessCost)
        mpc.hessCost = zeros(mpc.Nu+mpc.Nx);
    end
    
    [mpc.gradPerfQz,mpc.hessPerfTerm] = genLinOutGradHess(Qz,Cz,Dz,mpc.N,...
        mpc.N_ctr_hor,mpc.Nx,mpc.Nu,mpc.Nz,mpc.nx,mpc.nu,mpc.nz);

    mpc.hessCost = mpc.hessCost + mpc.hessPerfTerm;

end

% Linear Cost Term gradient Computation
if ~isempty(qz)

    mpc.gradPerfqz = genGenPerfLPGrad(qz,Cz,Dz,mpc.N,mpc.N_ctr_hor,...
                        mpc.Nx,mpc.Nu,mpc.nx,mpc.nu);

end
    
end