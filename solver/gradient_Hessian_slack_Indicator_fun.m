function [grad_fi_Ind,hess_fi_Ind] = gradient_Hessian_slack_Indicator_fun(cnstr,grad_fi_Ind,hess_fi_Ind)


if cnstr.min_limit

        grad_min_slack_Ind_x0 = grad_box_Ind(cnstr.min_slack_positivity_fi_x0,...
            cnstr.min_slack_positivity_grad);

        grad_fi_Ind = grad_fi_Ind + grad_min_slack_Ind_x0;

        hess_min_slack_Ind_x0 = hess_linear_Ind(cnstr.min_slack_positivity_fi_x0,...
            cnstr.min_slack_positivity_hess);

        hess_fi_Ind = hess_fi_Ind + hess_min_slack_Ind_x0;
end

if cnstr.max_limit
    
        grad_max_slack_Ind_x0 = grad_box_Ind(cnstr.max_slack_positivity_fi_x0,...
            cnstr.max_slack_positivity_grad);

        grad_fi_Ind = grad_fi_Ind + grad_max_slack_Ind_x0;

        hess_max_slack_Ind_x0 = hess_linear_Ind(cnstr.max_slack_positivity_fi_x0,...
                                                    cnstr.max_slack_positivity_hess);

        hess_fi_Ind = hess_fi_Ind + hess_max_slack_Ind_x0;
end

end