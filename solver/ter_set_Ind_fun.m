function [grad_ter,grad_ter_Ind_x0,hess_ter_Ind_x0] = ...
    ter_set_Ind_fun(mpc,x_ref)

% Gradient of Terminal Set Constraint
grad_ter = mpc.P2*(x_ref-mpc.s_ter); % grad_ter = -2*I*P*x_err -> negative I to be considered 
                                 % by substraction from global gradient 
if mpc.ter_constraint 
% fi_ter = e_xN'*P*e_xN - vi < 0
% grad_fi = 2*grad_e_xN*P*e_xN - [0;1]
% hess_fi = 2*grad_e_xN*P*grad_e_xN' (same as terminal cost function)

% 1. Indices for mapping
idx_xN = mpc.Nx+mpc.Nu-mpc.nx+1 : mpc.Nx+mpc.Nu;
idx_v = mpc.Nx+mpc.Nu+mpc.v_ter_global_index;
idx_active = [idx_xN, idx_v];
    
% Gradient of Indicator Function for Terminal Set Constraint
% -[0 -|grad_ter| 0 -|1|]/fi_ter_x0 negative cancels out
grad_ter_Ind_x0_sub = [grad_ter / mpc.fi_ter_x0;
                       1 / mpc.fi_ter_x0];

% Compute only the active sub-Hessian
hess_ter_Ind_x0_sub = grad_ter_Ind_x0_sub * grad_ter_Ind_x0_sub';
% Add the 2nd derivative of the constraint (-1/f * 2P)
hess_ter_Ind_x0_sub(1:mpc.nx, 1:mpc.nx) = ...
    hess_ter_Ind_x0_sub(1:mpc.nx, 1:mpc.nx) - mpc.P2 / mpc.fi_ter_x0;

% Map back to full vectors
% Allocate zeros and slot in ONLY the active blocks
grad_ter_Ind_x0 = zeros(mpc.Nx+mpc.Nu+mpc.Nv, 1);
grad_ter_Ind_x0(idx_active) = grad_ter_Ind_x0_sub;

hess_ter_Ind_x0 = zeros(mpc.Nx+mpc.Nu+mpc.Nv);
hess_ter_Ind_x0(idx_active, idx_active) = hess_ter_Ind_x0_sub;

else
    grad_ter_Ind_x0 = [];
    hess_ter_Ind_x0 = [];
end
 
end