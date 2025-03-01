% feasibility solver cost function is:
%
% f0 = v_feas + (x-x0)'Qfeas(x-x0)
%
% v_feas: feasibility slack variable
% (x-x0)'Qfeas(x-x0): quadritic term penalizing large deviations from the 
% original point
%
% for warm starting the mpc during initialization, the quadritic term
% should be ignored
%
% gradient computation: 
% Delta(v_feas) = [0,...,0,1]'
% Delta((x-x0)'Qfeas(x-x0)) = Qfeas*(x-x0)
%
% hessian computation:
% Hess((x-x0)'Qfeas(x-x0)) = Qfeas
function [grad_f0,hess_f0] = grad_hess_f0_feas(mpc,x,x0,n)

if mpc.warm_starting
    
    grad_f0 = [zeros(n,1);1];

    hess_f0 = [];

else

    grad_f0 = [zeros(n,1);1] + mpc.qfeas*[x(1:end-1)-x0;0];

    hess_f0 = diag([mpc.qfeas*ones(1,n), 0]);

end

end