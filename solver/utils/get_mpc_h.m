% Computes vectors: h = C*s + D*u Dsu*su + Dd*d
function [h_0,h,h_ter] = get_mpc_h( ...
    h_0,h,h_ter,s,u,su,dh,s_prev,u_prev,N, ...
    use_k0,use_s,use_u,use_su,use_d,use_ter, ...
    Ch_0,Dh_0,Dsuh_0,Ddh_0,Ch,Dh,Dsuh,Ddh,Ch_ter)

% k=0 (only if h depends on u)
if use_k0
    h_0(:) = Dh_0*u(:,1);
    if use_s
        h_0(:) = h_0(:) + Ch_0*s_prev;
    end
    if use_su
        h_0(:) = h_0(:) + Dsuh_0*u_prev;
    end
    if use_d
        h_0(:) = h_0(:) + Ddh_0*dh(:,1);
    end
end

h(:,:) = 0;
for k = 1:N-1

    if use_s
        h(:,k) = h(:,k) + Ch(:,:,k)*s(:,k);
    end
    if use_u
        h(:,k) = h(:,k) + Dh(:,:,k)*u(:,k+1);
    end
    if use_su
        h(:,k) = h(:,k) + Dsuh(:,:,k)*su(:,k);
    end
    if use_d
        h(:,k) = h(:,k) + Ddh(:,:,k)*dh(:,k+1);
    end
end

if use_ter
    h_ter(:) = Ch_ter*s(:,N);
end

end
