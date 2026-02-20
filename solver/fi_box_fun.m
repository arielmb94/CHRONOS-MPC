function [mpc,cnstr] = fi_box_fun(mpc,cnstr,x,N,n,v_feas)

if ~isempty(cnstr.min)

    % Let H be the number of time steps (horizon)
    H = N / n;

    % 1. Compute f_i for ALL time instances simultaneously
    cnstr.fi_min_x0 = repmat(cnstr.min, H, 1) - x;

    if cnstr.min_slack_nv > 0
        % Reshape and find the max violation across time steps for all constraints
        fi_mat = reshape(cnstr.fi_min_x0, n, H);
        max_fi_across_horizon = max(fi_mat, [], 2);

        % Compute candidate slacks for ALL variables at once
        candidate_v = max(mpc.slack_epsilon, max_fi_across_horizon + mpc.slack_epsilon);

        % Apply only to active constraints (multiplication masks out the inactive ones)
        cnstr.min_v = candidate_v .* cnstr.min_slack_active;

        % Map to global v vector
        mpc.v = mpc.v + cnstr.min_slack_map*cnstr.min_v;

        % 2. Apply the updated slack variables across all time instances
        cnstr.fi_min_x0 = cnstr.fi_min_x0 - repmat(cnstr.min_v, H, 1);
    end

    % for feasibility solver: fi - v_feas <= 0
    if v_feas
        cnstr.fi_min_x0 = cnstr.fi_min_x0 - v_feas;
    end

end

if ~isempty(cnstr.max)

    % Let H be the number of time steps (horizon)
    H = N / n;

    % 1. Compute f_i for ALL time instances simultaneously
    cnstr.fi_max_x0 = x - repmat(cnstr.max, H, 1);

    if cnstr.max_slack_nv > 0
        % Reshape and find the max violation across time steps for all constraints
        fi_mat = reshape(cnstr.fi_max_x0, n, H);
        max_fi_across_horizon = max(fi_mat, [], 2);

        % Compute candidate slacks for ALL variables at once
        candidate_v = max(mpc.slack_epsilon, max_fi_across_horizon + mpc.slack_epsilon);

        % Apply only to active constraints (multiplication masks out the inactive ones)
        cnstr.max_v = candidate_v .* cnstr.max_slack_active;

        % Map to global v vector
        mpc.v = mpc.v + cnstr.max_slack_map*cnstr.max_v;

        % 2. Apply the updated slack variables across all time instances
        cnstr.fi_max_x0 = cnstr.fi_max_x0 - repmat(cnstr.max_v, H, 1);
    end

    % for feasibility solver: fi - v_feas <= 0
    if v_feas
        cnstr.fi_max_x0 = cnstr.fi_max_x0 - v_feas;
    end

end

end