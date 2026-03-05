function [mpc,cnstr] = init_slack_max_condition(mpc,cnstr,slack_active_vector,...
                       N,n)

if length(slack_active_vector) == 1
    slack_active_vector = ones(n,1);
end

cnstr.max_slack_active = slack_active_vector;
nv = sum(slack_active_vector);
cnstr.max_slack_nv = nv;
% expand global slack vector
mpc.v = [mpc.v; zeros(nv,1)];
% create slack vector for max constraint
cnstr.max_v = zeros(n,1);
% map to get max_v to size nv
% vi = max_slack_local_map * max_v
cnstr.max_slack_local_map = create_shifted_identity(slack_active_vector);
% get index of vi in global v vector
cnstr.max_v_global_index = mpc.Nv+1:mpc.Nv+nv;

% increase constraint gradient vectors with slack variable contribution
grad_slack = -1 * create_shifted_identity(slack_active_vector);
cnstr.grad_max = [cnstr.grad_max;...
                  repmat(grad_slack,1,N/n)];

% slack variable gradient
max_slack_grad = [zeros(mpc.Nx+mpc.Nu+mpc.Nv,nv);...
                  eye(nv)];

% gradient/hess of constraint: -v<=0 (slack positivity constraint)
cnstr.max_slack_positivity_fi_x0 = zeros(nv,1);
cnstr.max_slack_positivity_grad = -1 * max_slack_grad;
[cnstr.max_slack_positivity_hess,mi] = genHessIneq(cnstr.max_slack_positivity_grad);
mpc.m = mpc.m+mi;

% Initialize Penalty term for new slack variables
if ~isempty(mpc.gradSlackqv)   
    % expand slack penalty term grad
    mpc.gradSlackqv = [mpc.gradSlackqv;
                       mpc.qv*ones(nv,1)];
else
    mpc.gradSlackqv = [zeros(mpc.Nx+mpc.Nu);
                        mpc.qv*ones(nv,1)];
end

% update global counter of slack variables
mpc.Nv = mpc.Nv+nv;

end