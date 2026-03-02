function [grad_fi_Ind,hess_fi_Ind] = gradient_Hessian_Indicator_fun(cnstr,grad_fi_Ind,hess_fi_Ind)


if ~isempty(cnstr.min)
    grad_min_Ind_x0 = grad_box_Ind(cnstr.fi_min_x0,cnstr.grad_min);

    grad_fi_Ind = grad_fi_Ind + grad_min_Ind_x0;

    hess_min_Ind_x0 = hess_linear_Ind(cnstr.fi_min_x0, cnstr.hess_min);

    hess_fi_Ind = hess_fi_Ind + hess_min_Ind_x0;

    if cnstr.min_slack_nv
        grad_min_slack_Ind_x0 = grad_box_Ind(cnstr.min_slack_positivity_fi_x0,...
            cnstr.min_slack_positivity_grad);

        grad_fi_Ind = grad_fi_Ind + grad_min_slack_Ind_x0;

        hess_min_slack_Ind_x0 = hess_linear_Ind(cnstr.min_slack_positivity_fi_x0,...
            cnstr.min_slack_positivity_hess);

        hess_fi_Ind = hess_fi_Ind + hess_min_slack_Ind_x0;
    end
end

if ~isempty(cnstr.max)
    grad_max_Ind_x0 = grad_box_Ind(cnstr.fi_max_x0,cnstr.grad_max);

    grad_fi_Ind = grad_fi_Ind + grad_max_Ind_x0;

    hess_max_Ind_x0 = hess_linear_Ind(cnstr.fi_max_x0,cnstr.hess_max);

    hess_fi_Ind = hess_fi_Ind + hess_max_Ind_x0;

    if cnstr.max_slack_nv
        grad_max_slack_Ind_x0 = grad_box_Ind(cnstr.max_slack_positivity_fi_x0,...
            cnstr.max_slack_positivity_grad);

        grad_fi_Ind = grad_fi_Ind + grad_max_slack_Ind_x0;

        hess_max_slack_Ind_x0 = hess_linear_Ind(cnstr.max_slack_positivity_fi_x0,...
                                                    cnstr.max_slack_positivity_hess);

        hess_fi_Ind = hess_fi_Ind + hess_max_slack_Ind_x0;
    end
end

end