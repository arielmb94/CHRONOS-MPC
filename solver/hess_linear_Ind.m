function hess = hess_linear_Ind(fi,hess_fi,...
                                active_set,active_indicator,n_i)

    % Hess(Phi) = sum((grad(fi)*grad(fi)^T)/fi^2) 
  
    m = length(hess_fi);
    n = size(hess_fi{1},1);

    hess = zeros(n);
    for i = 1:m
        
        if active_set
            index = mod(m,n_i);
            if ~index
                index = n_i;
            end
            if active_indicator(index)
                hess = hess + hess_fi{i}/(fi(i)^2);
            end
        else
            hess = hess + hess_fi{i}/(fi(i)^2);
        end
    end

end