function cnstr = fi_box_feas_fun(mpc,cnstr,x,N,n)

% Let H be the number of time steps (horizon)
H = N / n;

if ~isempty(cnstr.min)

    % 1. Compute f_i for ALL time instances simultaneously
    cnstr.fi_min_x0 = repmat(cnstr.min, H, 1) - x;

    if cnstr.min_slack_nv
        % 2. Apply slack variables across all time instances
        cnstr.fi_min_x0 = cnstr.fi_min_x0 - repmat(cnstr.min_v, H, 1);

        % update positivity constraints
        cnstr.min_slack_positivity_fi_x0 = mpc.slack_epsilon - ...
                                    cnstr.min_slack_local_map*cnstr.min_v;
    end

end

if ~isempty(cnstr.max)

    % 1. Compute f_i for ALL time instances simultaneously
    cnstr.fi_max_x0 = x - repmat(cnstr.max, H, 1);

    if cnstr.max_slack_nv
        % 2. Apply slack variables across all time instances
        cnstr.fi_max_x0 = cnstr.fi_max_x0 - repmat(cnstr.max_v, H, 1);

        % update positivity constraints
        cnstr.max_slack_positivity_fi_x0 = mpc.slack_epsilon - ...
                                    cnstr.max_slack_local_map*cnstr.max_v;
    end

end

end