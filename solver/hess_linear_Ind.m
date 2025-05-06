function hess = hess_linear_Ind(fi,hess_fi)

    % Hess(Phi) = sum((grad(fi)*grad(fi)^T)/fi^2) 
  
    m = length(hess_fi);
    n = size(hess_fi{1},1);

    hess = zeros(n);
    for i = 1:m

        hess = hess + hess_fi{i}/(fi(i)^2);

    end

end