%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Define constraint on the control variation between sampling instances:
%
%   du_min <= delta_u <= du_max
%
% You can define either lower or upper bounds alone, or both.
%
% Example uses:
%
%   - only lower bound: mpc = init_mpc_delta_u_cnstr(mpc,du_min,[])
%   - only upper bound: mpc = init_mpc_delta_u_cnstr(mpc,[],du_max)
%   - upper and lower bounds: mpc = init_mpc_delta_u_cnstr(mpc,du_min,du_max)
%
% In:
%   - mpc: CHRONOS mpc structure
%   - du_min (optional): nu column vector, lower bound constraint values on
%   the control variation between sampling instances
%   - du_max (optional): nu column vector, upper bound constraint values on
%   the control variation between sampling instances
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = init_mpc_delta_u_cnstr(mpc,du_min,du_max)

mpc.du_min = du_min;
mpc.du_max = du_max;

% Differential Control box constraints
[mpc.gradDeltaUmin,mpc.gradDeltaUmax] = genGradDeltaU(mpc.N_ctr_hor,...
                                                mpc.Nx,mpc.Nu,mpc.nx,mpc.nu);

if ~isempty(mpc.du_min)
    [mpc.hessDeltaUmin,mi] = genHessIneq(mpc.gradDeltaUmin);
    mpc.m = mpc.m+mi;
end
if ~isempty(mpc.du_max)
    [mpc.hessDeltaUmax,mi] = genHessIneq(mpc.gradDeltaUmax);
    mpc.m = mpc.m+mi;
end

end