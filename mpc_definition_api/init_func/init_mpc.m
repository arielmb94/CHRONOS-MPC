%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initializes CHRONOS mpc structure fields and solver hyperparmeters
%
% In:
%   - N: MPC prediction horizon
%   - N_ctr_hor (optional):  prediction horizon for control actions. If not
%   specified it is set equal to N
%
% Out:
%   - mpc: initialized CRHONOS mpc structure
%
% Example Use:
%
%   - Same control and prediction horizons:
%               mpc = init_mpc(N)
%   - Different control and prediction horizons:  
%               mpc = init_mpc(N,N_ctr_hor)
%
% Hyperparameters (can be modified manually after initialization of the mpc
% structure):
%
%   - mpc.t: interior-point method tradeoff parameter between cost function
%   minimization vs constraint satisfaction. Large t values give preference
%   to minimization of the cost function. Small values for t will make the
%   solver prefer feasibility and constraint safety.
%
%   - mpc.Beta: reduction step for each iteration of the feasibility line
%   search. Beta must be less than 1 and greater than 0. Values close to 1
%   ensure a smoother optimization solution between multiple mpc
%   iterations at the cost of increased line search iterations.
%
%   - mpc.min_l: if the line search step fall below min_l the following
%   iteration of the interior-point method will be cancelled. Allows to
%   quit the interior-point method quicker when the optimal solution is
%   close to the constraints limits.
%
%   - mpc.eps: interior-point method precision.
%
%   - mpc.max_iter: maximum allowed iterations of the interior-point method
%   solver
%
%   - mpc.t_feas: exactly as mpc.t, applied for the step 0 feasibility
%   solver. The step 0 solver allows to find a feasibile starting point for
%   the interior-point method when providded the initial guess lies
%   outside of the feasible region.
%
%   - mpc.qfeas: cost term to penalize large deviation on the solution of
%   the step 0 solver from the provided initial guess
%
%   - mpc.v0_feas: initial value for step 0 solver slack variable
%
%   - mpc.feas_lambda: multiplier in case step 0 starting slack variable
%   value is set too low. Must be larger than 1.
%
%   - mpc.max_feas_iter: maximum number of step 0 solver iterations 
%   allowed. If max_feas_iter is violated it is assumed the problem is 
%   unfeasible. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = init_mpc(N,N_ctr_hor)
arguments
    N = 2
    N_ctr_hor = 0
end

mpc.N = N;
if N_ctr_hor && N_ctr_hor > N
    mpc.N_ctr_hor = N;
elseif N_ctr_hor
    mpc.N_ctr_hor = N_ctr_hor;
else
    mpc.N_ctr_hor = N;
end
mpc.Qe = [];
mpc.Rdu = [];
mpc.Ru = [];
mpc.ru = [];
mpc.A = [];
mpc.B = [];
mpc.Bd = [];
mpc.C = [];
mpc.D = [];
mpc.Dd = [];
mpc.Qz = [];
mpc.qz = [];
mpc.Cz = [];
mpc.Dz = [];
mpc.Ddz = [];
mpc.Ch = [];
mpc.Dh = [];
mpc.Ddh = [];
mpc.nx = 0;
mpc.nu = 0;
mpc.nd = 0;
mpc.ny = 0;
mpc.ndz = 0;
mpc.nz = 0;
mpc.ndh = 0;
mpc.nh = 0;
mpc.Nx = 0;
mpc.Nu = 0;
mpc.Nd = 0;
mpc.Ny = 0;
mpc.Nz = 0;
mpc.Nh = 0;
mpc.Ndz = 0;
mpc.Ndh = 0;
mpc.x_min = [];
mpc.x_max = [];
mpc.x_ter_min = [];
mpc.x_ter_max = [];
mpc.u_min = [];
mpc.u_max = [];
mpc.du_min = [];
mpc.du_max = [];
mpc.y_min = [];
mpc.y_max = [];
mpc.h_min = [];
mpc.h_max = [];
mpc.Aeq = [];
mpc.beq = [];
mpc.hessCost = [];
mpc.gradErrQe = [];
mpc.hessErrTerm = [];
mpc.gradDiffCtlrRdu = [];
mpc.hessDiffCtrlTerm = [];
mpc.gradCtlrRu = [];
mpc.gradCtlrru = [];
mpc.hessCtrlTerm = [];
mpc.gradPerfQz = [];
mpc.gradPerfqz = [];
mpc.hessPerfTerm = [];
mpc.gradXmin = [];
mpc.gradXmax = [];
mpc.hessXmin = [];
mpc.hessXmax = [];
mpc.gradUmin = [];
mpc.gradUmax = [];
mpc.hessUmin = [];
mpc.hessUmax = [];
mpc.gradYmin = [];
mpc.gradYmax = [];
mpc.hessYmin = [];
mpc.hessYmax = [];
mpc.gradHmin = [];
mpc.gradHmax = [];
mpc.hessHmin = [];
mpc.hessHmax = [];
mpc.gradDeltaUmin = [];
mpc.gradDeltaUmax = [];
mpc.hessDeltaUmin = [];
mpc.hessDeltaUmax = [];
mpc.m = 0;
mpc.t = 50;
mpc.Beta = 0.75;
mpc.min_l = 0.99;
mpc.eps = 1e-4;
mpc.max_iter = 10;
mpc.ter_ingredients = 0;
mpc.ter_constraint = 0;
mpc.x_ref_is_y = 0;
mpc.P = [];
mpc.hessTerminalCost = [];
mpc.recompute_cost_hess = 0;
mpc.t_feas = 500;
mpc.v0_feas = 10;
mpc.qfeas = 1e-5;
mpc.warm_starting = 0;
mpc.feas_lambda = 100;
mpc.max_feas_iter = 10;
mpc.unfeasible = 0;
end