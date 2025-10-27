function [grad_ter,grad_ter_Ind_x0,hess_ter_Ind_x0] = ...
    ter_set_Ind_fun(x_ref,s_ter,fi_ter_x0,P2,Nx,Nu,nx,ter_constraint)

% Gradient of Terminal Set Constraint
grad_ter = P2*(x_ref-s_ter); % grad_ter = -2*I*P*x_err -> negative I to be considered 
                    % by substraction from global gradient 
if ter_constraint                    
% Gradient of Indicator Function for Terminal Set Constraint
grad_ter_Ind_x0 = zeros(Nx+Nu,1);
% -I * -grad_ter/fi_ter_x0 = grad_ter/fi_ter_x0
grad_ter_Ind_x0(end-(nx-1): end) = grad_ter/fi_ter_x0;

% Hessian of Indicator Function for Terminal Set Constraint
hess_ter_Ind_x0 = zeros(Nx+Nu,Nx+Nu);
hess_ter_Ind_x0(end-(nx-1): end, end-(nx-1): end) = ...
    (grad_ter*grad_ter')/(fi_ter_x0^2) - P2/fi_ter_x0;

else
    grad_ter_Ind_x0 = [];
    hess_ter_Ind_x0 = [];
end
 
end