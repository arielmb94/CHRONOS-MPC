function [mpc,cnstr] = init_slack_min_condition(mpc,cnstr,slack_active_vector,...
                       hard_limit_vector,N,n)

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

% compute vmax
if ~isempty(hard_limit_vector)
    cnstr.min_slack_hard_limit = 1;

    cnstr.min_slack_hard_limit_fi_x0 = zeros(nv,1);

    index = find(slack_active_vector);
    cnstr.min_slack_vmax = cnstr.min(index)-hard_limit_vector(index);

    % gradient/hess of constraint: v-vmax<=0 
    cnstr.min_slack_hard_limit_grad = cnstr.min_slack_index';
    [cnstr.min_slack_hard_limit_hess,mi] = genHessIneq(cnstr.min_slack_hard_limit_grad);
    mpc.m = mpc.m+mi;
else
    cnstr.min_slack_hard_limit = 0;
end  

% Initialize Penalty term for new slack variables
if ~isempty(mpc.gradSlackQv)   
    % expand slack penalty term grad
    mpc.gradSlackQv = [mpc.gradSlackQv;
                       zeros(nv,mpc.Nv)];
    % add columns for new slack
    mpc.gradSlackQv = [mpc.gradSlackQv min_slack_grad*mpc.Qv];

    % expand slack penalty term hess
    mpc.hessSlackTerm = [mpc.hessSlackTerm zeros(mpc.Nx+mpc.Nu+mpc.Nv,nv);
                         zeros(nv,mpc.Nx+mpc.Nu+mpc.Nv) eye(nv)*mpc.Qv];
else
    mpc.gradSlackQv = min_slack_grad*mpc.Qv;
    mpc.hessSlackTerm = [zeros(mpc.Nx+mpc.Nu) zeros(mpc.Nx+mpc.Nu,nv);
                         zeros(nv,mpc.Nx+mpc.Nu) eye(nv)*mpc.Qv];
end

% update global counter of slack variables
mpc.Nv = mpc.Nv+nv;

end