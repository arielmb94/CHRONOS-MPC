function mpc = f0_hess_reset_wrapper(mpc)
[mpc.recompute_cost_hess,...
    mpc.R_f0_0,mpc.Q_f0_k,mpc.R_f0_k,mpc.Y_f0_k,...
    mpc.Q_f0_ter] = f0_hess_reset_local(mpc.recompute_cost_hess,mpc.R_f0_0,...
                        mpc.Q_f0_k,mpc.R_f0_k,mpc.Y_f0_k,...
                        mpc.Q_f0_ter);
end

function [recompute_cost_hess,R_f0_0,Q_f0_k,R_f0_k,Y_f0_k,Q_f0_ter] = ...
    f0_hess_reset_local(recompute_cost_hess,R_f0_0,Q_f0_k,R_f0_k,Y_f0_k,Q_f0_ter)

recompute_cost_hess = false;
R_f0_0(:,:) = 0;
Q_f0_k(:,:,:) = 0;
R_f0_k(:,:,:) = 0;
Y_f0_k(:,:,:) = 0;
Q_f0_ter(:,:) = 0;
end
