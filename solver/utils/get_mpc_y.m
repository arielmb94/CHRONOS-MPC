% Computes output vectors of the form: y = Cs + Du + Dd
function [y_0,y,y_ter,err_0,err,err_ter] = get_mpc_y( ...
    y_0,y,y_ter,err_0,err,err_ter,s,u,d,r_0,r,r_ter,s_prev,N, ...
    tracking_cost,y_use_k0,y_use_s,y_use_u,y_use_d,y_use_ter, ...
    C_0,D_0,Dd_0,C,D,Dd,C_ter)

if y_use_k0
    y_0(:) = D_0*u(:,1);
    if y_use_s
        y_0(:) = y_0(:) + C_0*s_prev;
    end
    if y_use_d
        y_0(:) = y_0(:) + Dd_0*d(:,1);
    end

    if tracking_cost
        err_0(:) = r_0 - y_0;
    end
end

y(:,:) = 0;
for k = 1:N-1
    if y_use_s
        y(:,k) = y(:,k) + C(:,:,k)*s(:,k);
    end
    if y_use_u
        y(:,k) = y(:,k) + D(:,:,k)*u(:,k+1);
    end
    if y_use_d
        y(:,k) = y(:,k) + Dd(:,:,k)*d(:,k+1);
    end
end
if tracking_cost
    err(:,:) = r-y;
end

if y_use_ter
    y_ter(:) = C_ter*s(:,N);

    if tracking_cost
        err_ter(:) = r_ter - y_ter;
    end
end

end
