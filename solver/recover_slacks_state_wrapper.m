function mpc = recover_slacks_state_wrapper(mpc,s_cnstr)
[mpc.delta_g_k,mpc.delta_g_ter,mpc.delta_v_k,mpc.delta_v_ter] = ...
    recover_slacks_state_local(mpc.delta_g_k,mpc.delta_g_ter,mpc.delta_v_k,...
    mpc.delta_v_ter,mpc.N,mpc.delta_se,mpc.iS_ri_hat_k,mpc.iS_ri_hat_ter,...
    mpc.iS_k,mpc.iS_ter,mpc.g_k,mpc.g_ter,mpc.g2_k,mpc.g2_ter,...
    mpc.v2_k,mpc.v2_ter,mpc.rv_v2_k,mpc.rv_v2_ter,mpc.nx,...
    s_cnstr.min_limit,s_cnstr.max_limit,...
    s_cnstr.min_ineqRow_k,s_cnstr.max_ineqRow_k,...
    s_cnstr.min_row_v_k,s_cnstr.max_row_v_k,...
    s_cnstr.min_ineqRow_ter,s_cnstr.max_ineqRow_ter);
end

function [delta_g_k,delta_g_ter,delta_v_k,delta_v_ter] = ...
    recover_slacks_state_local(delta_g_k,delta_g_ter,delta_v_k,delta_v_ter,...
    N,delta_se,iS_ri_hat_k,iS_ri_hat_ter,iS_k,iS_ter,g_k,g_ter,g2_k,g2_ter,...
    v2_k,v2_ter,rv_v2_k,rv_v2_ter,nx,min_limit,max_limit,min_ineqRow_k,...
    max_ineqRow_k,min_row_v_k,max_row_v_k,min_ineqRow_ter,max_ineqRow_ter)
if ~isempty(min_limit) && ~isempty(max_limit)

    %k = 1,...N-1
    for k = 1:N-1
        for i = 1:nx
            row_min_i = min_ineqRow_k(i);
            row_max_i = max_ineqRow_k(i);

            mu_min_i = iS_ri_hat_k(row_min_i,k) - ...
                iS_k(row_min_i,k)*delta_se(i,k);
            mu_max_i = iS_ri_hat_k(row_max_i,k) + ...
                iS_k(row_max_i,k)*delta_se(i,k);

            g_min_i = g_k(row_min_i,k);
            g_max_i = g_k(row_max_i,k);
            g2_min_i = g2_k(row_min_i,k);
            g2_max_i = g2_k(row_max_i,k);

            delta_g_k(row_min_i,k) = g_min_i - g2_min_i*mu_min_i;
            delta_g_k(row_max_i,k) = g_max_i - g2_max_i*mu_max_i;

            row_v_min_i = min_row_v_k(i);
            row_v_max_i = max_row_v_k(i);
            v2_min_i = v2_k(row_v_min_i,k);
            v2_max_i = v2_k(row_v_max_i,k);

            delta_v_k(row_v_min_i,k) = ...
                -rv_v2_k(row_v_min_i,k) + v2_min_i*mu_min_i;
            delta_v_k(row_v_max_i,k) = ...
                -rv_v2_k(row_v_max_i,k) + v2_max_i*mu_max_i;
        end
    end

    %k = N
    for i = 1:nx
        row_min_i = min_ineqRow_ter(i);
        row_max_i = max_ineqRow_ter(i);

        mu_min_i = iS_ri_hat_ter(row_min_i) - ...
            iS_ter(row_min_i)*delta_se(i,N);
        mu_max_i = iS_ri_hat_ter(row_max_i) + ...
            iS_ter(row_max_i)*delta_se(i,N);

        g_min_i = g_ter(row_min_i);
        g_max_i = g_ter(row_max_i);
        g2_min_i = g2_ter(row_min_i);
        g2_max_i = g2_ter(row_max_i);

        delta_g_ter(row_min_i) = g_min_i - g2_min_i*mu_min_i;
        delta_g_ter(row_max_i) = g_max_i - g2_max_i*mu_max_i;

        v2_min_i = v2_ter(row_min_i);
        v2_max_i = v2_ter(row_max_i);

        delta_v_ter(row_min_i) = ...
            -rv_v2_ter(row_min_i) + v2_min_i*mu_min_i;
        delta_v_ter(row_max_i) = ...
            -rv_v2_ter(row_max_i) + v2_max_i*mu_max_i;
    end
elseif ~isempty(min_limit)

    %k = 1,...N-1
    for k = 1:N-1
        for i = 1:nx
            row_i = min_ineqRow_k(i);

            % A_i = -I <-- s
            mu_i = iS_ri_hat_k(row_i,k)-iS_k(row_i,k)*delta_se(i,k);

            g_i = g_k(row_i,k);
            g2_i = g2_k(row_i,k);

            delta_g_k(row_i,k) = g_i - g2_i*mu_i;

            row_v_i = min_row_v_k(i);
            v2_i = v2_k(row_v_i,k);

            delta_v_k(row_v_i,k) = -rv_v2_k(row_v_i,k) + v2_i*mu_i;
        end
    end

    %k = N
    for i = 1:nx
        row_i = min_ineqRow_ter(i);

        % A_i = -I <-- s
        mu_i = iS_ri_hat_ter(row_i)-iS_ter(row_i)*delta_se(i,N);

        g_i = g_ter(row_i);
        g2_i = g2_ter(row_i);
        delta_g_ter(row_i) = g_i - g2_i*mu_i;

        v2_i = v2_ter(row_i);
        delta_v_ter(row_i) = -rv_v2_ter(row_i) + v2_i*mu_i;
    end
elseif ~isempty(max_limit)

    %k = 1,...N-1
    for k = 1:N-1
        for i = 1:nx
            row_i = max_ineqRow_k(i);

            % A_i = I <-- s
            mu_i = iS_ri_hat_k(row_i,k)+iS_k(row_i,k)*delta_se(i,k);

            g_i = g_k(row_i,k);
            g2_i = g2_k(row_i,k);
            delta_g_k(row_i,k) = g_i - g2_i*mu_i;

            row_v_i = max_row_v_k(i);
            v2_i = v2_k(row_v_i,k);

            delta_v_k(row_v_i,k) = -rv_v2_k(row_v_i,k) + v2_i*mu_i;
        end
    end

    %k = N
    for i = 1:nx
        row_i = max_ineqRow_ter(i);

        % A_i = I <-- s
        mu_i = iS_ri_hat_ter(row_i)+iS_ter(row_i)*delta_se(i,N);

        g_i = g_ter(row_i);
        g2_i = g2_ter(row_i);
        delta_g_ter(row_i) = g_i - g2_i*mu_i;

        v2_i = v2_ter(row_i);
        delta_v_ter(row_i) = -rv_v2_ter(row_i) + v2_i*mu_i;
    end
end
end
