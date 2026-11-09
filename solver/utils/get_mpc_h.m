% Computes vectors: h = C*x + D*u Dxu*xu + Dd*d
function [h_0,h,h_ter] = get_mpc_h( ...
    h_0,h,h_ter,x,u,xu,dh,x_prev,u_prev,N, ...
    use_k0,use_x,use_u,use_xu,use_d,use_ter, ...
    Ch_0,Dh_0,Dxuh_0,Ddh_0,Ch,Dh,Dxuh,Ddh,Ch_ter)

% k=0 (only if h depends on u)
if use_k0
    h_0(:) = Dh_0*u(:,1);
    if use_x
        h_0(:) = h_0(:) + Ch_0*x_prev;
    end
    if use_xu
        h_0(:) = h_0(:) + Dxuh_0*u_prev;
    end
    if use_d
        h_0(:) = h_0(:) + Ddh_0*dh(:,1);
    end
end

h(:,:) = 0;
for k = 1:N-1

    if use_x
        h(:,k) = h(:,k) + Ch(:,:,k)*x(:,k);
    end
    if use_u
        h(:,k) = h(:,k) + Dh(:,:,k)*u(:,k+1);
    end
    if use_xu
        h(:,k) = h(:,k) + Dxuh(:,:,k)*xu(:,k);
    end
    if use_d
        h(:,k) = h(:,k) + Ddh(:,:,k)*dh(:,k+1);
    end
end

if use_ter
    h_ter(:) = Ch_ter*x(:,N);
end

end
