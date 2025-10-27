%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = update_mpc_sys_output(mpc,C,D,Dd,Qe,y_min,y_max)
%
% Allows to update all parameters related to the tracking output signal y 
% in a single function:
%
%   - Upadate output signal model: y = C * x + D * u + Dd * d
%   - Upadate weight Qe on the tracking error penalty term: 
%   (r - y)' * Qe * (r - y)
%   - Upadate contraint Limits on feedback signal: y_min <= y <= y_max
%
% Example uses:
%
%   - update only the output feedback signal model: 
%           mpc = update_mpc_sys_output(mpc,C,D,Dd)
%   - update only the input feedtrhough matrix of the output signal model: 
%           mpc = update_mpc_sys_output(mpc,[],D,[])
%   - update the feedback output signal model and constraint limits: 
%           mpc = update_mpc_sys_output(mpc,C,D,Dd,[],y_min,y_max)
%   - update only the weight on the tracking error penalty term: 
%           mpc = update_mpc_sys_output(mpc,[],[],[],Qe)
%
% In:
%   - mpc: CHRONOS mpc structure
%   - C (optional): ny x nx matrix, system output matrix
%   - D (optional): ny x nu matrix, input feedtrhough matrix.
%   - Dd (optional): ny x nd matrix, disturbance feedtrhough matrix.
%   - Qe (optional): ny x ny square matrix, weights for the quadratic
%   penalty on the tracking error
%   - y_min (optional): ny column vector, lower bound constraint values on 
%   the tracking signal
%   - y_max (optional): ny column vector, upper bound constraint values on 
%   the tracking signal
%
%   All arguments items which do not require to be updated can be passed as
%   an empty vector [].
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = update_mpc_sys_output(mpc,C,D,Dd,Qe,y_min,y_max)

update_gradients = 0;
update_cost_gradient = 0;

if ~isempty(C)
    mpc.C = C;
    update_gradients = 1;
end

if ~isempty(D)
    mpc.D = D;
    update_gradients = 1;
end

if ~isempty(Dd)   
    mpc.Dd = Dd;
end

if ~isempty(Qe)   
    mpc.Qe = Qe;
    update_cost_gradient = 1;
end

if ~isempty(mpc.y_cnstr) 
    if ~isempty(y_min)   
        mpc.y_cnstr.min = y_min;
    end
    
    if ~isempty(y_max)    
        mpc.y_cnstr.max = y_max;
    end
end

if update_cost_gradient || update_gradients

    % If tracking Cost exists, update gradients
    if ~isempty(mpc.Qe) || ~isempty(Qe)
        mpc = update_mpc_Tracking_cost(mpc,Qe);
    end
end

if update_gradients
    % If Outputs constraint exists, update box constraints gradients
    if ~isempty(mpc.y_cnstr) 
        [mpc.y_cnstr.grad_min(:,:),mpc.y_cnstr.grad_max(:,:)] = genGradY(mpc.C,mpc.D,mpc.N,mpc.N_ctr_hor,...
            mpc.Nx,mpc.Nu,mpc.Ny,mpc.nx,mpc.nu,mpc.ny);
    
        if ~isempty(mpc.y_cnstr.min)
            [mpc.y_cnstr.hess_min,~] = genHessIneq(mpc.y_cnstr.grad_min);
        end
        if ~isempty(mpc.y_cnstr.max)
            [mpc.y_cnstr.hess_max,~] = genHessIneq(mpc.y_cnstr.grad_max);
        end    
    end
end


end