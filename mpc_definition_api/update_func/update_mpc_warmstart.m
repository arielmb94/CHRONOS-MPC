%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   x0 = update_mpc_warmstart(x,mpc,fill_style)
%
% Provides alternative ways for warmstarting the optimization variables
% vector during runtime MPC executions.
% 
% The warmstarting style is set via the value of the variable fill_style.
% It controls the way in which the optimization variables vector is filled
% after warmstarting manipulations:
%
%   0: fills the warm-staterd optimization variables vector with 0s
%   1: fills the warm-staterd optimization variables vector by repeating 
%   the values from the prediction horizon step N-1.
%
% In:
%   - x: Nx + Nu column vector, optimization variables vector solution from
%   the last MPC iteration
%   - mpc: CHRONOS mpc structure
%   - fill_style (optional): warmstarting style variable
%
% Out:
%   - x0: warmstarted optimization variables vector for the follow-up MPC
%   iteration.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x0 = update_mpc_warmstart(x,mpc,fill_style)
arguments
    x
    mpc
    fill_style = 0;
end

    s_prev = zeros(mpc.nx,1); % s_prev is not used in warmstarting
    s = get_x(x,s_prev,mpc.nx,mpc.nu,mpc.N,mpc.N_ctr_hor,mpc.Nx);
    u = get_u(x,mpc.nx,mpc.nu,mpc.N_ctr_hor,mpc.Nu);

    % delete first state step and fill
    s_crop = s(mpc.nx+1:end);
    s_fill = fill_vec(s_crop,mpc.nx,mpc.Nx,fill_style);
    % delete first control step and fill
    u_crop = u(mpc.nu+1:end);
    u_fill = fill_vec(u_crop,mpc.nu,mpc.Nu,fill_style);

    x0 = arrange_opt_vec(s_fill,u_fill,mpc);

end

