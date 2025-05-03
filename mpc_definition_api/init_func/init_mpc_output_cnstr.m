%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = init_mpc_output_cnstr(mpc,y_min,y_max)
%
% Define constraints on the tracking signal y:
%
%   y_min <= y <= y_max
%
% In:
%   - mpc: CHRONOS mpc structure
%   - y_min (optional): ny column vector, lower bound constraint values on 
%   the tracking signal
%   - y_max (optional): ny column vector, upper bound constraint values on 
%   the tracking signal
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = init_mpc_output_cnstr(mpc,y_min,y_max)

mpc.y_min = y_min;
mpc.y_max = y_max;

% Outputs box constraints
[mpc.gradYmin,mpc.gradYmax] = genGradY(mpc.C,mpc.D,mpc.N,mpc.N_ctr_hor,...
    mpc.Nx,mpc.Nu,mpc.Ny,mpc.nx,mpc.nu,mpc.ny);

if ~isempty(mpc.y_min)
    [mpc.hessYmin,mi] = genHessIneq(mpc.gradYmin);
    mpc.m = mpc.m+mi;
end
if ~isempty(mpc.y_max)
    [mpc.hessYmax,mi] = genHessIneq(mpc.gradYmax);
    mpc.m = mpc.m+mi;
end

end