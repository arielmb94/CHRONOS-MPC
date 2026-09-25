function mpc = recover_slacks_output_wrapper(mpc,y_cnstr)
[mpc.delta_g_0,mpc.delta_g_k,mpc.delta_g_ter,mpc.delta_v_0,mpc.delta_v_k,...
    mpc.delta_v_ter] = recover_slacks_output_local(mpc.delta_g_0,mpc.delta_g_k,...
    mpc.delta_g_ter,mpc.delta_v_0,mpc.delta_v_k,mpc.delta_v_ter,mpc.N,...
    mpc.delta_u,mpc.delta_se,mpc.iS_ri_hat_0,mpc.iS_ri_hat_k,...
    mpc.iS_ri_hat_ter,mpc.iS_0,mpc.iS_k,mpc.iS_ter,mpc.g_0,mpc.g_k,...
    mpc.g_ter,mpc.g2_0,mpc.g2_k,mpc.g2_ter,mpc.v2_0,mpc.v2_k,mpc.v2_ter,...
    mpc.rv_v2_0,mpc.rv_v2_k,mpc.rv_v2_ter,mpc.D_0,mpc.C,mpc.D,mpc.C_ter,...
    mpc.s_col,y_cnstr.min_limit,y_cnstr.max_limit,y_cnstr.use_k0,...
    y_cnstr.use_s,y_cnstr.use_u,y_cnstr.use_ter,...
    y_cnstr.min_ineqRow_0,y_cnstr.max_ineqRow_0,...
    y_cnstr.min_row_v_0,y_cnstr.max_row_v_0,...
    y_cnstr.min_ineqRow_k,y_cnstr.max_ineqRow_k,...
    y_cnstr.min_row_v_k,y_cnstr.max_row_v_k,...
    y_cnstr.min_ineqRow_ter,y_cnstr.max_ineqRow_ter);
end

function [delta_g_0,delta_g_k,delta_g_ter,delta_v_0,delta_v_k,delta_v_ter] = ...
    recover_slacks_output_local(delta_g_0,delta_g_k,delta_g_ter,delta_v_0,...
    delta_v_k,delta_v_ter,N,delta_u,delta_se,iS_ri_hat_0,iS_ri_hat_k,...
    iS_ri_hat_ter,iS_0,iS_k,iS_ter,g_0,g_k,g_ter,g2_0,g2_k,g2_ter,v2_0,...
    v2_k,v2_ter,rv_v2_0,rv_v2_k,rv_v2_ter,D_0,C,D,C_ter,s_col,min_limit,...
    max_limit,use_k0,use_s,use_u,use_ter,min_ineqRow_0,max_ineqRow_0,...
    min_row_v_0,max_row_v_0,min_ineqRow_k,max_ineqRow_k,min_row_v_k,...
    max_row_v_k,min_ineqRow_ter,max_ineqRow_ter)
if ~isempty(min_limit) && ~isempty(max_limit)

    % k = 0
    if ~isempty(use_k0)
        delta_y_0 = D_0 * delta_u(:,1);

        mu_min = iS_ri_hat_0(min_ineqRow_0) - iS_0(min_ineqRow_0) .* delta_y_0;
        mu_max = iS_ri_hat_0(max_ineqRow_0) + iS_0(max_ineqRow_0) .* delta_y_0;

        delta_g_0(min_ineqRow_0) = g_0(min_ineqRow_0) - g2_0(min_ineqRow_0) .* mu_min;
        delta_g_0(max_ineqRow_0) = g_0(max_ineqRow_0) - g2_0(max_ineqRow_0) .* mu_max;

        delta_v_0(min_row_v_0) = -rv_v2_0(min_row_v_0) + v2_0(min_row_v_0) .* mu_min;
        delta_v_0(max_row_v_0) = -rv_v2_0(max_row_v_0) + v2_0(max_row_v_0) .* mu_max;
    end

    % k = 1,...,N-1
    for k = 1:N-1
        if ~isempty(use_s) && ~isempty(use_u)
            delta_y_k = C(:,:,k) * delta_se(s_col,k) + D(:,:,k) * delta_u(:,k+1);
        elseif ~isempty(use_s)
            delta_y_k = C(:,:,k) * delta_se(s_col,k);
        elseif ~isempty(use_u)
            delta_y_k = D(:,:,k) * delta_u(:,k+1);
        else
            delta_y_k = 0;
        end

        mu_min = iS_ri_hat_k(min_ineqRow_k,k) - iS_k(min_ineqRow_k,k) .* delta_y_k;
        mu_max = iS_ri_hat_k(max_ineqRow_k,k) + iS_k(max_ineqRow_k,k) .* delta_y_k;

        delta_g_k(min_ineqRow_k,k) = g_k(min_ineqRow_k,k) - g2_k(min_ineqRow_k,k) .* mu_min;
        delta_g_k(max_ineqRow_k,k) = g_k(max_ineqRow_k,k) - g2_k(max_ineqRow_k,k) .* mu_max;

        delta_v_k(min_row_v_k,k) = -rv_v2_k(min_row_v_k,k) + v2_k(min_row_v_k,k) .* mu_min;
        delta_v_k(max_row_v_k,k) = -rv_v2_k(max_row_v_k,k) + v2_k(max_row_v_k,k) .* mu_max;
    end

    % k = N
    if ~isempty(use_ter)
        delta_y_ter = C_ter * delta_se(s_col,N);

        mu_min = iS_ri_hat_ter(min_ineqRow_ter) - iS_ter(min_ineqRow_ter) .* delta_y_ter;
        mu_max = iS_ri_hat_ter(max_ineqRow_ter) + iS_ter(max_ineqRow_ter) .* delta_y_ter;

        delta_g_ter(min_ineqRow_ter) = g_ter(min_ineqRow_ter) - g2_ter(min_ineqRow_ter) .* mu_min;
        delta_g_ter(max_ineqRow_ter) = g_ter(max_ineqRow_ter) - g2_ter(max_ineqRow_ter) .* mu_max;

        delta_v_ter(min_ineqRow_ter) = -rv_v2_ter(min_ineqRow_ter) + v2_ter(min_ineqRow_ter) .* mu_min;
        delta_v_ter(max_ineqRow_ter) = -rv_v2_ter(max_ineqRow_ter) + v2_ter(max_ineqRow_ter) .* mu_max;
    end

elseif ~isempty(min_limit)

    % k = 0
    if ~isempty(use_k0)
        delta_y_0 = D_0 * delta_u(:,1);
        mu = iS_ri_hat_0(min_ineqRow_0) - iS_0(min_ineqRow_0) .* delta_y_0;

        delta_g_0(min_ineqRow_0) = g_0(min_ineqRow_0) - g2_0(min_ineqRow_0) .* mu;
        delta_v_0(min_row_v_0) = -rv_v2_0(min_row_v_0) + v2_0(min_row_v_0) .* mu;
    end

    % k = 1,...,N-1
    for k = 1:N-1
        if ~isempty(use_s) && ~isempty(use_u)
            delta_y_k = C(:,:,k) * delta_se(s_col,k) + D(:,:,k) * delta_u(:,k+1);
        elseif ~isempty(use_s)
            delta_y_k = C(:,:,k) * delta_se(s_col,k);
        elseif ~isempty(use_u)
            delta_y_k = D(:,:,k) * delta_u(:,k+1);
        else
            delta_y_k = 0;
        end

        mu = iS_ri_hat_k(min_ineqRow_k,k) - iS_k(min_ineqRow_k,k) .* delta_y_k;

        delta_g_k(min_ineqRow_k,k) = g_k(min_ineqRow_k,k) - g2_k(min_ineqRow_k,k) .* mu;
        delta_v_k(min_row_v_k,k) = -rv_v2_k(min_row_v_k,k) + v2_k(min_row_v_k,k) .* mu;
    end

    % k = N
    if ~isempty(use_ter)
        delta_y_ter = C_ter * delta_se(s_col,N);
        mu = iS_ri_hat_ter(min_ineqRow_ter) - iS_ter(min_ineqRow_ter) .* delta_y_ter;

        delta_g_ter(min_ineqRow_ter) = g_ter(min_ineqRow_ter) - g2_ter(min_ineqRow_ter) .* mu;
        delta_v_ter(min_ineqRow_ter) = -rv_v2_ter(min_ineqRow_ter) + v2_ter(min_ineqRow_ter) .* mu;
    end

elseif ~isempty(max_limit)

    % k = 0
    if ~isempty(use_k0)
        delta_y_0 = D_0 * delta_u(:,1);
        mu = iS_ri_hat_0(max_ineqRow_0) + iS_0(max_ineqRow_0) .* delta_y_0;

        delta_g_0(max_ineqRow_0) = g_0(max_ineqRow_0) - g2_0(max_ineqRow_0) .* mu;
        delta_v_0(max_row_v_0) = -rv_v2_0(max_row_v_0) + v2_0(max_row_v_0) .* mu;
    end

    % k = 1,...,N-1
    for k = 1:N-1
        if ~isempty(use_s) && ~isempty(use_u)
            delta_y_k = C(:,:,k) * delta_se(s_col,k) + D(:,:,k) * delta_u(:,k+1);
        elseif ~isempty(use_s)
            delta_y_k = C(:,:,k) * delta_se(s_col,k);
        elseif ~isempty(use_u)
            delta_y_k = D(:,:,k) * delta_u(:,k+1);
        else
            delta_y_k = 0;
        end

        mu = iS_ri_hat_k(max_ineqRow_k,k) + iS_k(max_ineqRow_k,k) .* delta_y_k;

        delta_g_k(max_ineqRow_k,k) = g_k(max_ineqRow_k,k) - g2_k(max_ineqRow_k,k) .* mu;
        delta_v_k(max_row_v_k,k) = -rv_v2_k(max_row_v_k,k) + v2_k(max_row_v_k,k) .* mu;
    end

    % k = N
    if ~isempty(use_ter)
        delta_y_ter = C_ter * delta_se(s_col,N);
        mu = iS_ri_hat_ter(max_ineqRow_ter) + iS_ter(max_ineqRow_ter) .* delta_y_ter;

        delta_g_ter(max_ineqRow_ter) = g_ter(max_ineqRow_ter) - g2_ter(max_ineqRow_ter) .* mu;
        delta_v_ter(max_ineqRow_ter) = -rv_v2_ter(max_ineqRow_ter) + v2_ter(max_ineqRow_ter) .* mu;
    end
end
end
