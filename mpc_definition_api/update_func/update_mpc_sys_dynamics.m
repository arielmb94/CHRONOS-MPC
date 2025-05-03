%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = update_mpc_sys_dynamics(mpc,A,B,Bd)
%
% Updates the matrices A, B and Bd of the internal MPC Discrete-Time system
% model:
% 
% x+ = A * x + B * u + Bd * d
% 
% and then recomputes the MPC equality constraints accordingly.
%
% Example uses:
%
%   - update only system matrix cost limits: 
%           mpc = update_mpc_sys_dynamics(mpc,A,[],[])
%   - update both system matrix and input matrox: 
%           mpc = update_mpc_sys_dynamics(mpc,[],B,[])
%   - update only input disturbance: 
%           mpc = update_mpc_sys_dynamics(mpc,[],[],Bd)
%
% In:
%   - mpc: CHRONOS mpc structure
%   - A (optional): nx x nx matrix, system matrix
%   - B (optional): nx x nu matrix, input matrix
%   - Bd (optional): nx x nd matrix, disturbance input matrix.
%
%   All arguments items which do not require to be updated can be passed as
%   an empty vector [].
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = update_mpc_sys_dynamics(mpc,A,B,Bd)

updateEqualities = 0;

if ~isempty(A)
    mpc.A = A;
    updateEqualities = 1;
end

if ~isempty(B)
    mpc.B = B;
    updateEqualities = 1;
end

if updateEqualities

    % A equality contraint 
    mpc.Aeq = genEqualities(mpc.A,mpc.B,mpc.N,mpc.N_ctr_hor,mpc.Nx,mpc.Nu,...
        mpc.nx,mpc.nu);
end

if ~isempty(Bd)   
    mpc.Bd = Bd;
end

end