function [mpc,cnstr] = init_slack_min_condition(mpc,cnstr,qv_min,N,n)

cnstr.has_slacks = 1;
% expand global slack vector
mpc.v = [mpc.v; zeros(n,1)];
% create slack vector for min constraint
cnstr.min_v = zeros(n,1);
% get index of vi in global v vector
cnstr.min_v_global_index = mpc.Nv+1:mpc.Nv+n;

% increase constraint gradient vectors with slack variable contribution
% fi-v<0 -> grad: grad_fi-eye(n)
grad_slack = -1 * eye(n);
cnstr.grad_min = [cnstr.grad_min;...
                  repmat(grad_slack,1,N/n)];

% slack variable gradient
min_slack_grad = [zeros(mpc.Nx+mpc.Nu+mpc.Nv,n);...
                  eye(n)];

% gradient/hess of constraint: -v<=0 (slack positivity constraint)
cnstr.min_slack_positivity_fi_x0 = zeros(n,1);
cnstr.min_slack_positivity_grad = -1 * min_slack_grad;
[cnstr.min_slack_positivity_hess,mi] = genHessIneq(cnstr.min_slack_positivity_grad);
mpc.m = mpc.m+mi;

% Initialize Penalty term for new slack variables
if isempty(qv_min)
    % if qv isnt defined, it is not initialized until build_chronos_mpc(), 
    % but we need to make space
    qv_min = mpc.qv*zeros(n,1);
elseif length(qv_min) == 1
    qv_min = qv_min*ones(n,1);
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
mpc.Nv = mpc.Nv+n;

end