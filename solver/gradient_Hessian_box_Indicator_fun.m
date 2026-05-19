function [grad_fi_Ind,hess_fi_Ind] = gradient_Hessian_box_Indicator_fun(cnstr,grad_fi_Ind,hess_fi_Ind)

if cnstr.min_limit

    grad_min_Ind_x0 = grad_box_Ind(cnstr.fi_min_x0,cnstr.grad_min);

    grad_fi_Ind = grad_fi_Ind + grad_min_Ind_x0;

    hess_min_Ind_x0 = hess_linear_Ind(cnstr.fi_min_x0, cnstr.hess_min);

    hess_fi_Ind = hess_fi_Ind + hess_min_Ind_x0;
end

if cnstr.max_limit

    grad_max_Ind_x0 = grad_box_Ind(cnstr.fi_max_x0,cnstr.grad_max);

    grad_fi_Ind = grad_fi_Ind + grad_max_Ind_x0;

    hess_max_Ind_x0 = hess_linear_Ind(cnstr.fi_max_x0,cnstr.hess_max);

    hess_fi_Ind = hess_fi_Ind + hess_max_Ind_x0;
end

end