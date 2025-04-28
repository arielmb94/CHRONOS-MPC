%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Define constraints on a custom user defined signal yi. The signal yi is
% defined as:
%
%   yi = Ci * x + Di * u + Ddi * di
%
% The user must define the signal yi by selecting appropiate values for the
% matrices Ci, Di, Ddi.
% 
% The constraints on yi are then of the form:
%
%   yi_min <= yi <= yi_max
%
% In:
%   - mpc: CHRONOS mpc structure
%   - yi_min: nyi column vector, lower bound constraint values on the user
%   defined signal yi
%   - yi_max: nyi column vector, upper bound constraint values on the user
%   defined signal yi
%   - Ci: nyi x nx matrix, states output matrix
%   - Di: nyi x nu matrix, input feedtrhough matrix
%   - Ddi: nyi x ndi matrix, disturbance feedtrhough matrix
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
% Example:
% We want to limit the variation of the MPC control action with respect a
% a given value u_star, e.g. we want to constraint yi = u - u_star.
% To achieve this we define yi by selecting:
%   - Ci = [0 0 ... 0]
%   - Di = [1]
%   - Ddi = [-1]
% which corresponds to yi = [0 0 ... 0] * x + [1] * u + [-1] * di
% The "disturbance" term on yi corresponds to u_star, to be introduced on 
% the appropiate field on mpc_solve() during runtime MPC execution.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = init_mpc_lin_custom_cnstr(mpc,yi_min,yi_max,Ci,Di,Ddi)

mpc.yi_min = yi_min;
mpc.yi_max = yi_max;

% General Inequality Matrix
mpc.Ci = Ci;
mpc.Di = Di;
mpc.Ddi = Ddi;

mpc.ndi = size(Ddi,2);  %number of disturbance inputs to general inequalities
mpc.nyi = size(Ci,1);  %number of general inequalities

if mpc.Di == 0
    mpc.Nyi = (mpc.N-1)*mpc.nyi;
else
    mpc.Nyi = mpc.N*mpc.nyi;
end

mpc.Ndi = mpc.N*mpc.ndi;

% General Inequalites box constraints
[mpc.gradYimin,mpc.gradYimax] = genGradY(mpc.Ci,mpc.Di,mpc.N,mpc.N_ctr_hor,...
                                mpc.Nx,mpc.Nu,mpc.Nyi,mpc.nx,mpc.nu,mpc.nyi);

if ~isempty(mpc.yi_min)
    [mpc.hessYimin,mi] = genHessIneq(mpc.gradYimin);
    mpc.m = mpc.m+mi;
end
if ~isempty(mpc.yi_max)
    [mpc.hessYimax,mi] = genHessIneq(mpc.gradYimax);
    mpc.m = mpc.m+mi;
end

end