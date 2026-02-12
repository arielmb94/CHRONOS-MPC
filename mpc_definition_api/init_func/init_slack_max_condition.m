function [mpc,cnstr] = init_slack_max_condition(mpc,cnstr,slack_active_vector,...
                       hard_limit_vector,N,n)

if length(slack_active_vector) == 1
    slack_active_vector = ones(n,1);
end

cnstr.max_slack_map = create_shifted_identity(slack_active_vector)';

nv = sum(slack_active_vector);
cnstr.max_slack_nv = nv;

cnstr.max_slack_positivity_fi_x0 = zeros(nv,1);

grad_slack = -1 * create_shifted_identity(slack_active_vector);

% increase gradient vectors with slack variable gradient
cnstr.grad_max = [cnstr.grad_max;...
                  repmat(grad_slack,1,N/n)];

% slack variable index in optimization vector
cnstr.max_slack_index = [zeros(nv,mpc.Nx+mpc.Nu+mpc.Nv) eye(nv)];

% gradient/hess of constraint: -v<=0 
cnstr.max_slack_positivity_grad = -1 * cnstr.max_slack_index';
[cnstr.max_slack_positivity_hess,mi] = genHessIneq(cnstr.max_slack_positivity_grad);
mpc.m = mpc.m+mi;

% compute vmax
if ~isempty(hard_limit_vector)
    cnstr.max_slack_hard_limit = 1;
    index = find(slack_active_vector);
    cnstr.max_slack_vmax = hard_limit_vector(index)-cnstr.max(index);

    cnstr.max_slack_hard_limit_fi_x0 = zeros(nv,1);

    % gradient/hess of constraint: v-vmax<=0 
    cnstr.max_slack_hard_limit_grad = cnstr.max_slack_index';
    [cnstr.max_slack_hard_limit_hess,mi] = genHessIneq(cnstr.max_slack_hard_limit_grad);
    mpc.m = mpc.m+mi;
else
    cnstr.max_slack_hard_limit = 0;
end

% Initialize Penalty term for new slack variables
if ~isempty(mpc.gradSlackQv)   
    % expand slack penalty term grad
    mpc.gradSlackQv = [mpc.gradSlackQv;
                       zeros(nv,mpc.Nv)];
    % add columns for new slack
    mpc.gradSlackQv = [mpc.gradSlackQv cnstr.max_slack_index'*mpc.Qv];

    % expand slack penalty term hess
    mpc.hessSlackTerm = [mpc.hessSlackTerm zeros(mpc.Nx+mpc.Nu+mpc.Nv,nv);
                         zeros(nv,mpc.Nx+mpc.Nu+mpc.Nv) eye(nv)*mpc.Qv];
else
    mpc.gradSlackQv = cnstr.max_slack_index'*mpc.Qv;
    mpc.hessSlackTerm = [zeros(mpc.Nx+mpc.Nu) zeros(mpc.Nx+mpc.Nu,nv);
                         zeros(nv,mpc.Nx+mpc.Nu) eye(nv)*mpc.Qv];
end

% update global counter of slack variables
mpc.Nv = mpc.Nv+nv;

end