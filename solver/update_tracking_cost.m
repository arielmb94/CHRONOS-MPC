function mpc = update_tracking_cost(mpc,update_tracking)

[mpc.update_tracking,mpc.recompute_cost_hess,mpc.grad_u_E_0,mpc.grad_s_E, ...
 mpc.grad_u_E,mpc.grad_s_E_ter,mpc.R_E_0,mpc.Q_E,mpc.R_E,mpc.Y_E,mpc.Q_E_ter] = ...
 update_tracking_cost_local(update_tracking,mpc.recompute_cost_hess, ...
 mpc.grad_u_E_0,mpc.grad_s_E,mpc.grad_u_E,mpc.grad_s_E_ter,mpc.R_E_0, ...
 mpc.Q_E,mpc.R_E,mpc.Y_E,mpc.Q_E_ter,mpc.y_use_k0,mpc.y_use_s, ...
 mpc.y_use_u,mpc.y_use_ter,mpc.Qe_0,mpc.Qe,mpc.Qe_ter,mpc.D_0,mpc.C, ...
 mpc.D,mpc.C_ter,mpc.N);
end

function [update_tracking,recompute_cost_hess,grad_u_E_0,grad_s_E,grad_u_E, ...
 grad_s_E_ter,R_E_0,Q_E,R_E,Y_E,Q_E_ter] = update_tracking_cost_local( ...
 update_tracking,recompute_cost_hess,grad_u_E_0,grad_s_E,grad_u_E,grad_s_E_ter, ...
 R_E_0,Q_E,R_E,Y_E,Q_E_ter,y_use_k0,y_use_s,y_use_u,y_use_ter,Qe_0,Qe, ...
 Qe_ter,D_0,C,D,C_ter,N)

update_tracking = false;
recompute_cost_hess = true;
if ~isempty(y_use_k0)
 grad_u_E_0(:,:) = D_0'*Qe_0;
 R_E_0(:,:) = grad_u_E_0*D_0;
end
if ~isempty(y_use_s) && ~isempty(y_use_u)
 for k=1:N-1
  grad_s_E(:,:,k) = C(:,:,k)'*Qe(:,:,k);
  grad_u_E(:,:,k) = D(:,:,k)'*Qe(:,:,k);
  Q_E(:,:,k) = grad_s_E(:,:,k)*C(:,:,k);
  R_E(:,:,k) = grad_u_E(:,:,k)*D(:,:,k);
  Y_E(:,:,k) = grad_u_E(:,:,k)*C(:,:,k);
 end
elseif ~isempty(y_use_s)
 for k=1:N-1
  grad_s_E(:,:,k) = C(:,:,k)'*Qe(:,:,k);
  Q_E(:,:,k) = grad_s_E(:,:,k)*C(:,:,k);
 end
elseif ~isempty(y_use_u)
 for k=1:N-1
  grad_u_E(:,:,k) = D(:,:,k)'*Qe(:,:,k);
  R_E(:,:,k) = grad_u_E(:,:,k)*D(:,:,k);
 end
end
if ~isempty(y_use_ter)
 grad_s_E_ter(:,:) = C_ter'*Qe_ter;
 Q_E_ter(:,:) = grad_s_E_ter*C_ter;
end
end
