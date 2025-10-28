% n_i is size of
function grad_Ind_x0 = grad_box_Ind(fi,grad_fi,...
                                    active_set,active_indicator,n_i)

    [n,m] = size(grad_fi);

    grad_Ind_x0 = zeros(n,1);
    for i = 1:m

        if active_set
            index = mod(m,n_i);
            if ~index
                index = n_i;
            end
            if active_indicator(index)
                grad_Ind_x0 = grad_Ind_x0 - grad_fi(:,i)/fi(i);
            end
        else
            grad_Ind_x0 = grad_Ind_x0 - grad_fi(:,i)/fi(i);
        end
    end

end