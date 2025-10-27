%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = update_mpc_Control_cost(mpc,Ru,ru)
%
% Modifies the weights Ru and ru for the quadratic and linear control 
% penalty terms on the control action. The function then updates the MPC 
% gradients and Hessians accordingly.
%
% In:
%   - mpc: CHRONOS mpc structure
%   - Ru (optional): nu x nu square matrix, weights for the quadratic
%   penalty term on the control action.
%   - ru (optional): nu column vector, weights for the linear penalty term
%   on the control action. IMPORTANT: Use linear penalties only in the case
%   that the control action takes positive values only.
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = update_mpc_Control_cost(mpc,Ru,ru)

if ~isempty(Ru)
    
    mpc.Ru = Ru;
    
    [mpc.gradCtlrRu(:,:),mpc.hessCtrlTerm(:,:)] = genControlGradHess(Ru,mpc.N_ctr_hor,...
        mpc.Nx,mpc.Nu,mpc.nx,mpc.nu);
    
    mpc.recompute_cost_hess = 1;

end

if ~isempty(ru)

    mpc.ru = ru;

    mpc.gradCtlrru(:,:) = genControlLPGrad(ru,mpc.N_ctr_hor,mpc.Nx,mpc.Nu,...
                    mpc.nx,mpc.nu);

end

end