function [grad_ter_Ind_x0,hess_ter_Ind_x0] = ...
    ter_set_feas_Ind_fun(x_ref,s_ter,fi_ter_x0,P,nx,n)
% optimization variables are stored as [u0 x1 u1 x2 ... xn un xn+1 v]^T
% v is feasibility problem slack variable

% Gradient of Terminal Set Constraint
grad_ter = [2*P*(x_ref-s_ter);1]; % grad_ter = [-I*P*x_err;-1] = -[P*x_err;1] 

% Gradient of Indicator Function for Terminal Set Constraint
grad_ter_Ind_x0 = zeros(n,1);

% Negative of Indicator gradient cancels with negative of terminal ser
% gradient
% - -[P*x_err;1]/fi_ter_x0 = grad_ter/fi_ter_x0
grad_ter_Ind_x0(end-(nx) : end) = grad_ter/fi_ter_x0;

% Hessian of Indicator Function for Terminal Set Constraint
hess_ter_Ind_x0 = zeros(n);

Hess_fi_ter = zeros(nx+1);
Hess_fi_ter(1:nx,1:nx) = 2*P;

hess_ter_Ind_x0(end-(nx) : end, end-(nx) : end) = ...
    (grad_ter*grad_ter')/(fi_ter_x0^2) - Hess_fi_ter/fi_ter_x0;

end