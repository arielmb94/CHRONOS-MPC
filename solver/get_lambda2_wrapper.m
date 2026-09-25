function lambda2 = get_lambda2_wrapper(mpc)
lambda2 = get_lambda2_local(mpc.delta_u,mpc.delta_se,...
    mpc.delta_g_0,mpc.delta_g_k,mpc.delta_g_ter,...
    mpc.delta_v_0,mpc.delta_v_k,mpc.delta_v_ter,...
    mpc.g_0,mpc.g_k,mpc.g_ter,mpc.v_0,mpc.v_k,mpc.v_ter,...
    mpc.grad_u_f0_0,mpc.grad_se_f0_k,mpc.grad_u_f0_k,mpc.grad_se_f0_ter,...
    mpc.t,mpc.grad_qv_0,mpc.grad_qv_k,mpc.grad_qv_ter,mpc.N);
end

function lambda2 = get_lambda2_local(...
    delta_u,delta_se,delta_g_0,delta_g_k,delta_g_ter,...
    delta_v_0,delta_v_k,delta_v_ter,g_0,g_k,g_ter,v_0,v_k,v_ter,...
    grad_u_f0_0,grad_se_f0_k,grad_u_f0_k,grad_se_f0_ter,...
    t,grad_qv_0,grad_qv_k,grad_qv_ter,N)
lambda2 = 0;

% k = 0: f0, g, v
lambda2 = lambda2-grad_u_f0_0'*delta_u(:,1);
if ~isempty(g_0)
    lambda2 = lambda2+(1./g_0)'*delta_g_0;
end
if ~isempty(v_0)
    lambda2 = lambda2-(t*grad_qv_0-1./v_0)'*delta_v_0;
end

% k = 1:N-1: f0, g, v
for k = 1:N-1
    lambda2 = lambda2-grad_se_f0_k(:,k)'*delta_se(:,k);
    lambda2 = lambda2-grad_u_f0_k(:,k)'*delta_u(:,k+1);
    if ~isempty(g_k)
        lambda2 = lambda2+(1./g_k(:,k))'*delta_g_k(:,k);
    end
    if ~isempty(v_k)
        lambda2 = lambda2-(t*grad_qv_k(:,k)-1./v_k(:,k))'*delta_v_k(:,k);
    end
end

% k = N: f0, g, v
lambda2 = lambda2-grad_se_f0_ter'*delta_se(:,N);
if ~isempty(g_ter)
    lambda2 = lambda2+(1./g_ter)'*delta_g_ter;
end
if ~isempty(v_ter)
    lambda2 = lambda2-(t*grad_qv_ter-1./v_ter)'*delta_v_ter;
end
end
