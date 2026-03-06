function [mpc,cnstr] = init_slack_min_condition(mpc,cnstr,slack_active_vector,...
                       qv_min,N,n)

if length(slack_active_vector) == 1
    slack_active_vector = ones(n,1);
end

cnstr.min_slack_active = slack_active_vector;
nv = sum(slack_active_vector);
cnstr.min_slack_nv = nv;
% expand global slack vector
mpc.v = [mpc.v; zeros(nv,1)];
% create slack vector for min constraint
cnstr.min_v = zeros(n,1);
% map to get min_v to size nv
% vi = min_slack_local_map * min_v
cnstr.min_slack_local_map = create_shifted_identity(slack_active_vector);
% get index of vi in global v vector
cnstr.min_v_global_index = mpc.Nv+1:mpc.Nv+nv;

% increase constraint gradient vectors with slack variable contribution
grad_slack = -1 * create_shifted_identity(slack_active_vector);
cnstr.grad_min = [cnstr.grad_min;...
                  repmat(grad_slack,1,N/n)];

% slack variable gradient
min_slack_grad = [zeros(mpc.Nx+mpc.Nu+mpc.Nv,nv);...
                  eye(nv)];

% gradient/hess of constraint: -v<=0 (slack positivity constraint)
cnstr.min_slack_positivity_fi_x0 = zeros(nv,1);
cnstr.min_slack_positivity_grad = -1 * min_slack_grad;
[cnstr.min_slack_positivity_hess,mi] = genHessIneq(cnstr.min_slack_positivity_grad);
mpc.m = mpc.m+mi;

% Initialize Penalty term for new slack variables
if isempty(qv_min)
    qv_min = mpc.qv*ones(nv,1);
elseif length(qv_min) == 1
    qv_min = qv_min*ones(nv,1);
else
    qv_min = cnstr.min_slack_local_map*qv_min;
end

if ~isempty(mpc.gradSlackqv)   
    % expand slack penalty term grad
    mpc.gradSlackqv = [mpc.gradSlackqv;
                       qv_min];
else
    mpc.gradSlackqv = [zeros(mpc.Nx+mpc.Nu,1);
                        qv_min];
end

% update global counter of slack variables
mpc.Nv = mpc.Nv+nv;

end