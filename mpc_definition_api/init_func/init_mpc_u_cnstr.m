%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = init_mpc_u_cnstr(mpc,u_min,u_max)
%
% Define limits on the control action:
%
%   u_min <= u <= u_max
%
% You can define either lower or upper bounds alone, or both.
%
% Example uses:
%
%   - only lower bound: mpc = init_mpc_u_cnstr(mpc,u_min,[])
%   - only upper bound: mpc = init_mpc_u_cnstr(mpc,[],u_max)
%   - upper and lower bounds: mpc = init_mpc_u_cnstr(mpc,u_min,u_max)
%
% In:
%   - mpc: CHRONOS mpc structure
%   - u_min (optional): nu column vector, lower bound constraint values on
%   the control action
%   - u_max (optional): nu column vector, upper bound constraint values on
%   the control action
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = init_mpc_u_cnstr(mpc,u_min,u_max)

mpc.u_min = u_min;
mpc.u_max = u_max;

% Control box constraints
[mpc.gradUmin,mpc.gradUmax] = genGradU(mpc.N_ctr_hor,...
                                    mpc.Nx,mpc.Nu,mpc.nx,mpc.nu);

if ~isempty(mpc.u_min)
    [mpc.hessUmin,mi] = genHessIneq(mpc.gradUmin);
    mpc.m = mpc.m+mi;
end
if ~isempty(mpc.u_max)
    [mpc.hessUmax,mi] = genHessIneq(mpc.gradUmax);
    mpc.m = mpc.m+mi;
end

end