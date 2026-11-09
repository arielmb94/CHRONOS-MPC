function mpc = mpc_h_wrapper(mpc,h_cnstr,x_prev,u_prev)
[mpc.h_0,mpc.h,mpc.h_ter] = mpc_h_local(mpc.h_0,mpc.h,mpc.h_ter,mpc.x,mpc.u, ...
    mpc.xu,mpc.dh,x_prev,u_prev,mpc.N,h_cnstr.use_k0,h_cnstr.use_x, ...
    h_cnstr.use_u,h_cnstr.use_xu,h_cnstr.use_d,h_cnstr.use_ter, ...
    mpc.Ch_0,mpc.Dh_0,mpc.Dxuh_0,mpc.Ddh_0,mpc.Ch,mpc.Dh,mpc.Dxuh,mpc.Ddh,mpc.Ch_ter);
end

function [h_0,h,h_ter] = mpc_h_local(h_0,h,h_ter,x,u,xu,dh,x_prev,u_prev,N, ...
    use_k0,use_x,use_u,use_xu,use_d,use_ter,Ch_0,Dh_0,Dxuh_0,Ddh_0,Ch,Dh,Dxuh,Ddh,Ch_ter)
if ~isempty(use_k0)
    h_0(:) = Dh_0*u(:,1);
    if ~isempty(use_x), h_0(:) = h_0(:) + Ch_0*x_prev; end
    if ~isempty(use_xu), h_0(:) = h_0(:) + Dxuh_0*u_prev; end
    if ~isempty(use_d), h_0(:) = h_0(:) + Ddh_0*dh(:,1); end
end
h(:,:) = 0;
for k = 1:N-1
    if ~isempty(use_x), h(:,k) = h(:,k) + Ch(:,:,k)*x(:,k); end
    if ~isempty(use_u), h(:,k) = h(:,k) + Dh(:,:,k)*u(:,k+1); end
    if ~isempty(use_xu), h(:,k) = h(:,k) + Dxuh(:,:,k)*xu(:,k); end
    if ~isempty(use_d), h(:,k) = h(:,k) + Ddh(:,:,k)*dh(:,k+1); end
end
if ~isempty(use_ter), h_ter(:) = Ch_ter*x(:,N); end
end
