% Updates fi (fi<0) for box constraints with slack -> fi-vi<0
% Computes positivity constraint fi for slacks epsilon<v -> fi := epsilon-v 
function cnstr = fi_box_slack_compute(mpc,cnstr,N,n)

% Let H be the number of time steps (horizon)
H = N / n;

if cnstr.min_limit
    % Apply slack variables across all time instances
    cnstr.fi_min_x0 = cnstr.fi_min_x0 - repmat(cnstr.min_v, H, 1);

    % update positivity constraints
    cnstr.min_slack_positivity_fi_x0 = mpc.slack_epsilon - cnstr.min_v;
end

if cnstr.max_limit
    % Apply slack variables across all time instances
    cnstr.fi_max_x0 = cnstr.fi_max_x0 - repmat(cnstr.max_v, H, 1);

    % update positivity constraints
    cnstr.max_slack_positivity_fi_x0 = mpc.slack_epsilon - cnstr.max_v;
end

end