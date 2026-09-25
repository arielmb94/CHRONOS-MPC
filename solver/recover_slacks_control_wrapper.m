function mpc = recover_slacks_control_wrapper(mpc,u_cnstr)
[mpc.delta_g_0,mpc.delta_g_k] = recover_slacks_control_local(...
    mpc.delta_g_0,mpc.delta_g_k,mpc.N,mpc.delta_u,mpc.iS_ri_hat_0,...
    mpc.iS_ri_hat_k,mpc.g_0,mpc.g_k,mpc.g2_0,mpc.g2_k,mpc.nu,...
    u_cnstr.min_limit,u_cnstr.max_limit,...
    u_cnstr.min_ineqRow_0,u_cnstr.max_ineqRow_0,...
    u_cnstr.min_ineqRow_k,u_cnstr.max_ineqRow_k);
end

function [delta_g_0,delta_g_k] = recover_slacks_control_local(delta_g_0,delta_g_k,...
    N,delta_u,iS_ri_hat_0,iS_ri_hat_k,g_0,g_k,g2_0,g2_k,nu,min_limit,...
    max_limit,min_ineqRow_0,max_ineqRow_0,min_ineqRow_k,max_ineqRow_k)
if ~isempty(min_limit) && ~isempty(max_limit)
    %k = 0
    for i = 1:nu
        row_min_i = min_ineqRow_0(i);
        row_max_i = max_ineqRow_0(i);

        g_min_i = g_0(row_min_i);
        g_max_i = g_0(row_max_i);
        g2_min_i = g2_0(row_min_i);
        g2_max_i = g2_0(row_max_i);

        delta_g_0(row_min_i) = g_min_i - ...
            g2_min_i*iS_ri_hat_0(row_min_i) + delta_u(i,1);
        delta_g_0(row_max_i) = g_max_i - ...
            g2_max_i*iS_ri_hat_0(row_max_i) - delta_u(i,1);
    end

    %k = 1,...N-1
    for k = 1:N-1
        for i = 1:nu
            row_min_i = min_ineqRow_k(i);
            row_max_i = max_ineqRow_k(i);

            g_min_i = g_k(row_min_i,k);
            g_max_i = g_k(row_max_i,k);
            g2_min_i = g2_k(row_min_i,k);
            g2_max_i = g2_k(row_max_i,k);

            delta_g_k(row_min_i,k) = g_min_i - ...
                g2_min_i*iS_ri_hat_k(row_min_i,k) + delta_u(i,k+1);
            delta_g_k(row_max_i,k) = g_max_i - ...
                g2_max_i*iS_ri_hat_k(row_max_i,k) - delta_u(i,k+1);
        end
    end
elseif ~isempty(min_limit)
    %k = 0
    for i = 1:nu
        row_i = min_ineqRow_0(i);

        % A_i = -I <-- u
        g_i = g_0(row_i);
        g2_i = g2_0(row_i);
        delta_g_0(row_i) = g_i - g2_i*iS_ri_hat_0(row_i) + delta_u(i,1);
    end

    %k = 1,...N-1
    for k = 1:N-1
        for i = 1:nu
            row_i = min_ineqRow_k(i);

            % A_i = -I <-- u
            g_i = g_k(row_i,k);
            g2_i = g2_k(row_i,k);
            delta_g_k(row_i,k) = g_i - ...
                g2_i*iS_ri_hat_k(row_i,k) + delta_u(i,k+1);
        end
    end
elseif ~isempty(max_limit)
    %k = 0
    for i = 1:nu
        row_i = max_ineqRow_0(i);

        % A_i = I <-- u
        g_i = g_0(row_i);
        g2_i = g2_0(row_i);
        delta_g_0(row_i) = g_i - g2_i*iS_ri_hat_0(row_i) - delta_u(i,1);
    end
    %k = 1,...N-1
    for k = 1:N-1
        for i = 1:nu
            row_i = max_ineqRow_k(i);

            % A_i = I <-- u
            g_i = g_k(row_i,k);
            g2_i = g2_k(row_i,k);
            delta_g_k(row_i,k) = g_i - ...
                g2_i*iS_ri_hat_k(row_i,k) - delta_u(i,k+1);
        end
    end
end
end
