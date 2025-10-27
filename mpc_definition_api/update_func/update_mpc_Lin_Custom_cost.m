%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = update_mpc_Lin_Custom_cost(mpc,Cz,Dz,Ddz,Qz,qz)
%
% Updates the parameters for the user-defined cost term and re-computes the
% MPC gradients and Hessians accordingly.
%
% The function can be used to update the model of the user-defined signal 
% z:
%
%   z = Cz * x + Dz * u + Ddz * dz
%
% by updating the matrices Cz, Dz and Ddz.
%
% The function can also be called to update the weight values for the
% quadratic Ru and linear ru penalty terms of the MPC cost functions:
%
%   J += z'*Qz*z + qz*z
%
% Example uses:
%
%   - update only quadratic cost penalty: 
%           mpc = update_mpc_Lin_Custom_cost(mpc,[],[],[],Qz)
%   - update only user-defined signal z model: 
%           mpc = update_mpc_Lin_Custom_cost(mpc,Cz,Dz,Ddz)
%   - update only input feedthrough Dz matrix and linear cost: 
%           mpc = update_mpc_Lin_Custom_cost(mpc,[],Di,[],[],qz)
%
% In:
%   - mpc: CHRONOS mpc structure
%   - Cz (optional): nz x nx matrix, state output matrix
%   - Dz (optional): nz x nu matrix, input feedtrhough output matrix
%   - Ddz (optional): nz x ndz matrix, disturbance feedtrhough output 
%   matrix
%   - Qz (optional): nz x nz square matrix, weights for the quadratic
%   penalty term on the user defined signal z
%   - qz (optional): nz column vector, weights for the linear penalty term
%   on the user defined signal z.
%
%   All arguments items which do not require to be updated can be passed as
%   an empty vector [].
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = update_mpc_Lin_Custom_cost(mpc,Cz,Dz,Ddz,Qz,qz)

updatelp = 0;
updateqp = 0;

if ~isempty(Cz)
    mpc.Cz = Cz;

    updatelp = 1;
    updateqp = 1;
end

if ~isempty(Dz)
    mpc.Dz = Dz;

    updatelp = 1;
    updateqp = 1;
end

if ~isempty(Ddz)   
    mpc.Ddz = Ddz;
end

if ~isempty(Qz)   
    mpc.Qz = Qz;

    updateqp = 1;
end

if ~isempty(qz)   
    mpc.qz = qz;

    updatelp = 1;
end

% Update Quatric cost term gradient and Hessian
if ~isempty(mpc.Qz) && updateqp

    mpc.recompute_cost_hess = 1;

    [mpc.gradPerfQz(:,:),mpc.hessPerfTerm(:,:)] = genLinOutGradHess(mpc.Qz, ...
        mpc.Cz,mpc.Dz,mpc.N,mpc.N_ctr_hor,mpc.Nx,mpc.Nu,mpc.Nz,...
        mpc.nx,mpc.nu,mpc.nz);

end

% Update Linear cost term gradient
if ~isempty(mpc.qz) && updatelp

    mpc.gradPerfqz(:,:) = genGenPerfLPGrad(mpc.qz,mpc.Cz,mpc.Dz,...
        mpc.N,mpc.N_ctr_hor,mpc.Nx,mpc.Nu,mpc.nx,mpc.nu);

end
    
end