function grad_Ind_x0 = grad_box_Ind(fi,grad_fi)

    [n,m] = size(grad_fi);

    grad_Ind_x0 = zeros(n,1);
    for i = 1:m

        grad_Ind_x0 = grad_Ind_x0 - grad_fi(:,i)/fi(i);

    end

end