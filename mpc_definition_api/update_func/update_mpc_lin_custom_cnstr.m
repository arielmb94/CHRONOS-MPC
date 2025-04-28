%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Updates the parameters for the user-defined custom constraint and
% re-computes the MPC gradients and Hessians accordingly.
% 
% The function can be used to update the model of the user-defined signal 
% yi:
%
%   yi = Ci * x + Di * u + Ddi * di
%
% by updating the matrices Ci, Di and Ddi.
%
% The function can also be called to update the constraint limits:
%
%   yi_min <= yi <= yi_max
%
% Example uses:
%
%   - update only constraint limits: 
%           mpc = init_mpc_delta_u_cnstr(mpc,[],[],[],yi_min,yi_max)
%   - update only user-defined signal yi model: 
%           mpc = init_mpc_delta_u_cnstr(mpc,Ci,Di,Ddi)
%   - update only input feedthrough Di matrix : 
%           mpc = init_mpc_delta_u_cnstr(mpc,[],Di,[])
%
% In:
%   - mpc: CHRONOS mpc structure
%   - Ci (optional): nyi x nx matrix, state output matrix
%   - Di (optional): nyi x nu matrix, input feedtrhough matrix
%   - Ddi (optional): nyi x ndi matrix, disturbance feedtrhough matrix
%   - yi_min (optional): nyi column vector, lower bound constraint values
%   on the user defined signal yi
%   - yi_max (optional): nyi column vector, upper bound constraint values 
%   on the user defined signal yi
%
%   All arguments items which do not require to be updated can be passed as
%   an empty vector [].
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = update_mpc_lin_custom_cnstr(mpc,Ci,Di,Ddi,...
                yi_min,yi_max)
arguments
    mpc
    Ci
    Di
    Ddi
    yi_min = []
    yi_max = []
end

update_grads = 0;

if ~isempty(Ci)
    mpc.Ci = Ci;
    update_grads = 1;
end

if ~isempty(Di)
    mpc.Di = Di;
    update_grads = 1;
end

if ~isempty(Ddi)   
    mpc.Ddi = Ddi;
end

if ~isempty(yi_min)   
    mpc.yi_min = yi_min;
end

if ~isempty(yi_max)    
    mpc.yi_max = yi_max;
end

% General Inequalites box constraints
if update_grads
    
    [mpc.gradYimin,mpc.gradYimax] = genGradY(mpc.Ci,mpc.Di,mpc.N,mpc.N_ctr_hor,...
        mpc.Nx,mpc.Nu,mpc.Nyi,mpc.nx,mpc.nu,mpc.nyi);
    
    if ~isempty(mpc.yi_min)
        [mpc.hessYimin,~] = genHessIneq(mpc.gradYimin);
    end
    if ~isempty(mpc.yi_max)
        [mpc.hessYimax,~] = genHessIneq(mpc.gradYimax);
    end

end

end