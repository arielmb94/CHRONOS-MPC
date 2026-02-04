function [mpc,cnstr] = init_slack_max_condition(mpc,cnstr,slack_active_vector,...
    N,n)

if length(slack_active_vector) == 1
    slack_active_vector = ones(n,1);
end

nv = sum(slack_active_vector);
cnstr.max_slack_nv = nv;

grad_slack = -1 * create_shifted_identity(slack_active_vector);

% increase gradient vectors with slack variable gradient
cnstr.grad_max = [cnstr.grad_max;...
                  repmat(grad_slack,1,N/n)];

% slack variable index in optimization vector
cnstr.max_slack_index = [zeros(nv,mpc.Nx+mpc.Nu+mpc.Nv) eye(nv)];

% update global counter of slack variables
mpc.Nv = mpc.Nv+nv;

end