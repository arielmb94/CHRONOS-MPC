function mpc = recover_slacks_custom_wrapper(mpc,h_cnstr)
[mpc.delta_g_0,mpc.delta_g_k,mpc.delta_g_ter,mpc.delta_v_0,mpc.delta_v_k,...
    mpc.delta_v_ter] = recover_slacks_custom_local(mpc.delta_g_0,mpc.delta_g_k,...
    mpc.delta_g_ter,mpc.delta_v_0,mpc.delta_v_k,mpc.delta_v_ter,mpc.N,...
    mpc.delta_u,mpc.delta_se,mpc.iS_ri_hat_0,mpc.iS_ri_hat_k,...
    mpc.iS_ri_hat_ter,mpc.iS_0,mpc.iS_k,mpc.iS_ter,mpc.g_0,mpc.g_k,...
    mpc.g_ter,mpc.g2_0,mpc.g2_k,mpc.g2_ter,mpc.v2_0,mpc.v2_k,mpc.v2_ter,...
    mpc.rv_v2_0,mpc.rv_v2_k,mpc.rv_v2_ter,mpc.Dh_0,mpc.Ch,mpc.Dsuh,mpc.Dh,...
    mpc.Ch_ter,mpc.s_col,mpc.su_col,h_cnstr.min_limit,...
    h_cnstr.max_limit,h_cnstr.use_k0,h_cnstr.use_s,...
    h_cnstr.use_su,h_cnstr.use_u,h_cnstr.use_ter,...
    h_cnstr.min_ineqRow_0,h_cnstr.max_ineqRow_0,...
    h_cnstr.min_row_v_0,h_cnstr.max_row_v_0,...
    h_cnstr.min_ineqRow_k,h_cnstr.max_ineqRow_k,...
    h_cnstr.min_row_v_k,h_cnstr.max_row_v_k,...
    h_cnstr.min_ineqRow_ter,h_cnstr.max_ineqRow_ter);
end

function [delta_g_0,delta_g_k,delta_g_ter,delta_v_0,delta_v_k,delta_v_ter] = ...
    recover_slacks_custom_local(delta_g_0,delta_g_k,delta_g_ter,delta_v_0,...
    delta_v_k,delta_v_ter,N,delta_u,delta_se,iS_ri_hat_0,iS_ri_hat_k,...
    iS_ri_hat_ter,iS_0,iS_k,iS_ter,g_0,g_k,g_ter,g2_0,g2_k,g2_ter,v2_0,...
    v2_k,v2_ter,rv_v2_0,rv_v2_k,rv_v2_ter,Dh_0,Ch,Dsuh,Dh,Ch_ter,s_col,...
    su_col,min_limit,max_limit,use_k0,use_s,use_su,use_u,use_ter,...
    min_ineqRow_0,max_ineqRow_0,min_row_v_0,max_row_v_0,min_ineqRow_k,...
    max_ineqRow_k,min_row_v_k,max_row_v_k,min_ineqRow_ter,max_ineqRow_ter)

if ~isempty(min_limit) && ~isempty(max_limit)

    % k = 0
    if ~isempty(use_k0)
        delta_h_0 = Dh_0 * delta_u(:,1);

        mu_min = iS_ri_hat_0(min_ineqRow_0) - iS_0(min_ineqRow_0) .* delta_h_0;
        mu_max = iS_ri_hat_0(max_ineqRow_0) + iS_0(max_ineqRow_0) .* delta_h_0;

        delta_g_0(min_ineqRow_0) = g_0(min_ineqRow_0) - g2_0(min_ineqRow_0) .* mu_min;
        delta_g_0(max_ineqRow_0) = g_0(max_ineqRow_0) - g2_0(max_ineqRow_0) .* mu_max;

        delta_v_0(min_row_v_0) = -rv_v2_0(min_row_v_0) + v2_0(min_row_v_0) .* mu_min;
        delta_v_0(max_row_v_0) = -rv_v2_0(max_row_v_0) + v2_0(max_row_v_0) .* mu_max;
    end

    % k = 1,...,N-1
    for k = 1:N-1
        if ~isempty(use_s) && ~isempty(use_su) && ~isempty(use_u)
            delta_h_k = Ch(:,:,k) * delta_se(s_col,k) + ...
                Dsuh(:,:,k) * delta_se(su_col,k) + Dh(:,:,k) * delta_u(:,k+1);
        elseif ~isempty(use_s) && ~isempty(use_su)
            delta_h_k = Ch(:,:,k) * delta_se(s_col,k) + Dsuh(:,:,k) * delta_se(su_col,k);
        elseif ~isempty(use_s) && ~isempty(use_u)
            delta_h_k = Ch(:,:,k) * delta_se(s_col,k) + Dh(:,:,k) * delta_u(:,k+1);
        elseif ~isempty(use_su) && ~isempty(use_u)
            delta_h_k = Dsuh(:,:,k) * delta_se(su_col,k) + Dh(:,:,k) * delta_u(:,k+1);
        elseif ~isempty(use_s)
            delta_h_k = Ch(:,:,k) * delta_se(s_col,k);
        elseif ~isempty(use_su)
            delta_h_k = Dsuh(:,:,k) * delta_se(su_col,k);
        elseif ~isempty(use_u)
            delta_h_k = Dh(:,:,k) * delta_u(:,k+1);
        end

        mu_min = iS_ri_hat_k(min_ineqRow_k,k) - iS_k(min_ineqRow_k,k) .* delta_h_k;
        mu_max = iS_ri_hat_k(max_ineqRow_k,k) + iS_k(max_ineqRow_k,k) .* delta_h_k;

        delta_g_k(min_ineqRow_k,k) = g_k(min_ineqRow_k,k) - g2_k(min_ineqRow_k,k) .* mu_min;
        delta_g_k(max_ineqRow_k,k) = g_k(max_ineqRow_k,k) - g2_k(max_ineqRow_k,k) .* mu_max;

        delta_v_k(min_row_v_k,k) = -rv_v2_k(min_row_v_k,k) + v2_k(min_row_v_k,k) .* mu_min;
        delta_v_k(max_row_v_k,k) = -rv_v2_k(max_row_v_k,k) + v2_k(max_row_v_k,k) .* mu_max;
    end

    % k = N
    if ~isempty(use_ter)
        delta_h_ter = Ch_ter * delta_se(s_col,N);

        mu_min = iS_ri_hat_ter(min_ineqRow_ter) - iS_ter(min_ineqRow_ter) .* delta_h_ter;
        mu_max = iS_ri_hat_ter(max_ineqRow_ter) + iS_ter(max_ineqRow_ter) .* delta_h_ter;

        delta_g_ter(min_ineqRow_ter) = g_ter(min_ineqRow_ter) - g2_ter(min_ineqRow_ter) .* mu_min;
        delta_g_ter(max_ineqRow_ter) = g_ter(max_ineqRow_ter) - g2_ter(max_ineqRow_ter) .* mu_max;

        delta_v_ter(min_ineqRow_ter) = -rv_v2_ter(min_ineqRow_ter) + v2_ter(min_ineqRow_ter) .* mu_min;
        delta_v_ter(max_ineqRow_ter) = -rv_v2_ter(max_ineqRow_ter) + v2_ter(max_ineqRow_ter) .* mu_max;
    end

elseif ~isempty(min_limit)

    % k = 0
    if ~isempty(use_k0)
        delta_h_0 = Dh_0 * delta_u(:,1);
        mu = iS_ri_hat_0(min_ineqRow_0) - iS_0(min_ineqRow_0) .* delta_h_0;

        delta_g_0(min_ineqRow_0) = g_0(min_ineqRow_0) - g2_0(min_ineqRow_0) .* mu;
        delta_v_0(min_row_v_0) = -rv_v2_0(min_row_v_0) + v2_0(min_row_v_0) .* mu;
    end

    % k = 1,...,N-1
    for k = 1:N-1
        if ~isempty(use_s) && ~isempty(use_su) && ~isempty(use_u)
            delta_h_k = Ch(:,:,k) * delta_se(s_col,k) + ...
                Dsuh(:,:,k) * delta_se(su_col,k) + Dh(:,:,k) * delta_u(:,k+1);
        elseif ~isempty(use_s) && ~isempty(use_su)
            delta_h_k = Ch(:,:,k) * delta_se(s_col,k) + Dsuh(:,:,k) * delta_se(su_col,k);
        elseif ~isempty(use_s) && ~isempty(use_u)
            delta_h_k = Ch(:,:,k) * delta_se(s_col,k) + Dh(:,:,k) * delta_u(:,k+1);
        elseif ~isempty(use_su) && ~isempty(use_u)
            delta_h_k = Dsuh(:,:,k) * delta_se(su_col,k) + Dh(:,:,k) * delta_u(:,k+1);
        elseif ~isempty(use_s)
            delta_h_k = Ch(:,:,k) * delta_se(s_col,k);
        elseif ~isempty(use_su)
            delta_h_k = Dsuh(:,:,k) * delta_se(su_col,k);
        elseif ~isempty(use_u)
            delta_h_k = Dh(:,:,k) * delta_u(:,k+1);
        end

        mu = iS_ri_hat_k(min_ineqRow_k,k) - iS_k(min_ineqRow_k,k) .* delta_h_k;

        delta_g_k(min_ineqRow_k,k) = g_k(min_ineqRow_k,k) - g2_k(min_ineqRow_k,k) .* mu;
        delta_v_k(min_row_v_k,k) = -rv_v2_k(min_row_v_k,k) + v2_k(min_row_v_k,k) .* mu;
    end

    % k = N
    if ~isempty(use_ter)
        delta_h_ter = Ch_ter * delta_se(s_col,N);
        mu = iS_ri_hat_ter(min_ineqRow_ter) - iS_ter(min_ineqRow_ter) .* delta_h_ter;

        delta_g_ter(min_ineqRow_ter) = g_ter(min_ineqRow_ter) - g2_ter(min_ineqRow_ter) .* mu;
        delta_v_ter(min_ineqRow_ter) = -rv_v2_ter(min_ineqRow_ter) + v2_ter(min_ineqRow_ter) .* mu;
    end

elseif ~isempty(max_limit)

    % k = 0
    if ~isempty(use_k0)
        delta_h_0 = Dh_0 * delta_u(:,1);
        mu = iS_ri_hat_0(max_ineqRow_0) + iS_0(max_ineqRow_0) .* delta_h_0;

        delta_g_0(max_ineqRow_0) = g_0(max_ineqRow_0) - g2_0(max_ineqRow_0) .* mu;
        delta_v_0(max_row_v_0) = -rv_v2_0(max_row_v_0) + v2_0(max_row_v_0) .* mu;
    end

    % k = 1,...,N-1
    for k = 1:N-1
        if ~isempty(use_s) && ~isempty(use_su) && ~isempty(use_u)
            delta_h_k = Ch(:,:,k) * delta_se(s_col,k) + ...
                Dsuh(:,:,k) * delta_se(su_col,k) + Dh(:,:,k) * delta_u(:,k+1);
        elseif ~isempty(use_s) && ~isempty(use_su)
            delta_h_k = Ch(:,:,k) * delta_se(s_col,k) + Dsuh(:,:,k) * delta_se(su_col,k);
        elseif ~isempty(use_s) && ~isempty(use_u)
            delta_h_k = Ch(:,:,k) * delta_se(s_col,k) + Dh(:,:,k) * delta_u(:,k+1);
        elseif ~isempty(use_su) && ~isempty(use_u)
            delta_h_k = Dsuh(:,:,k) * delta_se(su_col,k) + Dh(:,:,k) * delta_u(:,k+1);
        elseif ~isempty(use_s)
            delta_h_k = Ch(:,:,k) * delta_se(s_col,k);
        elseif ~isempty(use_su)
            delta_h_k = Dsuh(:,:,k) * delta_se(su_col,k);
        elseif ~isempty(use_u)
            delta_h_k = Dh(:,:,k) * delta_u(:,k+1);
        end

        mu = iS_ri_hat_k(max_ineqRow_k,k) + iS_k(max_ineqRow_k,k) .* delta_h_k;

        delta_g_k(max_ineqRow_k,k) = g_k(max_ineqRow_k,k) - g2_k(max_ineqRow_k,k) .* mu;
        delta_v_k(max_row_v_k,k) = -rv_v2_k(max_row_v_k,k) + v2_k(max_row_v_k,k) .* mu;
    end

    % k = N
    if ~isempty(use_ter)
        delta_h_ter = Ch_ter * delta_se(s_col,N);
        mu = iS_ri_hat_ter(max_ineqRow_ter) + iS_ter(max_ineqRow_ter) .* delta_h_ter;

        delta_g_ter(max_ineqRow_ter) = g_ter(max_ineqRow_ter) - g2_ter(max_ineqRow_ter) .* mu;
        delta_v_ter(max_ineqRow_ter) = -rv_v2_ter(max_ineqRow_ter) + v2_ter(max_ineqRow_ter) .* mu;
    end

end
end
