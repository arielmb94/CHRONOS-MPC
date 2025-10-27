function hessIneq = updateHessIneq(grad_fi,hessIneq)

    % Hess(Phi) = sum((grad(fi)*grad(fi)^T)/fi^2)   
    m = size(grad_fi,2);

    for i = 1:m
        hessIneq{i} = grad_fi(:,i)*grad_fi(:,i)';
    end

end