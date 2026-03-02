function [cnstr,v] = fi_box_slack_update(mpc,cnstr,x,N,n)

% Let H be the number of time steps (horizon)
H = N / n;

v = zeros(mpc.Nv,1);
v(:) = mpc.v;

if ~isempty(cnstr.min)

    % 1. Compute f_i for ALL time instances simultaneously
    cnstr.fi_min_x0 = repmat(cnstr.min, H, 1) - x;

    if cnstr.min_slack_nv
        % Reshape and find the max violation across time steps for all constraints
        fi_mat = reshape(cnstr.fi_min_x0, n, H);
        max_fi_across_horizon = max(fi_mat, [], 2);

        % Compute candidate slacks for ALL variables at once
        candidate_v = max([mpc.slack_epsilon*ones(n,1),...
                           max_fi_across_horizon + mpc.slack_epsilon,...
                           cnstr.min_v],[], 2);

        % Apply only to active constraints (multiplication masks out the inactive ones)
        cnstr.min_v = candidate_v .* cnstr.min_slack_active;

        % 2. Apply slack variables across all time instances
        cnstr.fi_min_x0 = cnstr.fi_min_x0 - repmat(cnstr.min_v, H, 1);

        % Map to global v vector
        vi = cnstr.min_slack_local_map*cnstr.min_v;
        v(cnstr.min_v_global_index) = vi;

        % update positivity constraints
        cnstr.min_slack_positivity_fi_x0 = mpc.slack_epsilon - vi;
    end

end

if ~isempty(cnstr.max)

    % 1. Compute f_i for ALL time instances simultaneously
    cnstr.fi_max_x0 = x - repmat(cnstr.max, H, 1);

    if cnstr.max_slack_nv
        % Reshape and find the max violation across time steps for all constraints
        fi_mat = reshape(cnstr.fi_max_x0, n, H);
        max_fi_across_horizon = max(fi_mat, [], 2);

        % Compute candidate slacks for ALL variables at once
        candidate_v = max([mpc.slack_epsilon*ones(n,1),...
                           max_fi_across_horizon + mpc.slack_epsilon,...
                           cnstr.max_v],[], 2);

        % Apply only to active constraints (multiplication masks out the inactive ones)
        cnstr.max_v = candidate_v .* cnstr.max_slack_active;

        % 2. Apply slack variables across all time instances
        cnstr.fi_max_x0 = cnstr.fi_max_x0 - repmat(cnstr.max_v, H, 1);

        % Map to global v vector
        vi = cnstr.max_slack_local_map*cnstr.max_v;
        v(cnstr.max_v_global_index) = vi;

        % update positivity constraints
        cnstr.max_slack_positivity_fi_x0 = mpc.slack_epsilon - vi;
    end

end

end