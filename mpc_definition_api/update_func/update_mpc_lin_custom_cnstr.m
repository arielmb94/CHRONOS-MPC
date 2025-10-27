%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    mpc = update_mpc_lin_custom_cnstr(mpc,Ch,Dh,Ddh,h_min,h_max)
%
% Updates the parameters for the user-defined custom constraint and
% re-computes the MPC gradients and Hessians accordingly.
% 
% The function can be used to update the model of the user-defined signal 
% h:
%
%   h = Ch * x + Dh * u + Ddh * dh
%
% by updating the matrices Ch, Dh and Ddh.
%
% The function can also be called to update the constraint limits:
%
%   h_min <= h <= h_max
%
% Example uses:
%
%   - update only constraint limits: 
%           mpc = init_mpc_delta_u_cnstr(mpc,[],[],[],h_min,h_max)
%   - update only user-defined signal h model: 
%           mpc = init_mpc_delta_u_cnstr(mpc,Ch,Dh,Ddh)
%   - update only input feedthrough Dh matrix : 
%           mpc = init_mpc_delta_u_cnstr(mpc,[],Dh,[])
%
% In:
%   - mpc: CHRONOS mpc structure
%   - Ch (optional): nh x nx matrix, state output matrix
%   - Dh (optional): nh x nu matrix, input feedtrhough matrix
%   - Ddh (optional): nh x ndh matrix, disturbance feedtrhough matrix
%   - h_min (optional): nh column vector, lower bound constraint values
%   on the user defined signal h
%   - h_max (optional): nh column vector, upper bound constraint values 
%   on the user defined signal h
%
%   All arguments items which do not require to be updated can be passed as
%   an empty vector [].
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = update_mpc_lin_custom_cnstr(mpc,Ch,Dh,Ddh,h_min,h_max)

update_grads = 0;

if ~isempty(Ch)
    mpc.Ch = Ch;
    update_grads = 1;
end

if ~isempty(Dh)
    mpc.Dh = Dh;
    update_grads = 1;
end

if ~isempty(Ddh)   
    mpc.Ddh = Ddh;
end

if ~isempty(h_min)   
    mpc.h_cnstr.min = h_min;
end

if ~isempty(h_max)    
    mpc.h_cnstr.max = h_max;
end

% General Inequalites box constraints
if update_grads
    
    [mpc.h_cnstr.grad_min(:,:),mpc.h_cnstr.grad_max(:,:)] = genGradY(mpc.Ch,mpc.Dh,mpc.N,mpc.N_ctr_hor,...
        mpc.Nx,mpc.Nu,mpc.Nh,mpc.nx,mpc.nu,mpc.nh);
    
    if ~isempty(mpc.h_cnstr.min)
        [mpc.h_cnstr.hess_min,~] = genHessIneq(mpc.h_cnstr.grad_min);
    end
    if ~isempty(mpc.h_cnstr.max)
        [mpc.h_cnstr.hess_max,~] = genHessIneq(mpc.h_cnstr.grad_max);
    end

end

end