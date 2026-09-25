function mpc = reduced_kkt_control_rate_wrapper(mpc,du_cnstr)
[mpc.ru_0,mpc.rse_k,mpc.ru_k,mpc.R_0,mpc.Q_k,mpc.R_k,mpc.Y_k] = reduced_kkt_control_rate_local(mpc.N, mpc.nu, mpc.R_0, mpc.Q_k, mpc.R_k, mpc.Y_k, mpc.ru_0, mpc.rse_k, mpc.ru_k, du_cnstr.min_limit, du_cnstr.max_limit, du_cnstr.min_ineqRow_0, du_cnstr.max_ineqRow_0, du_cnstr.min_ineqRow_k, du_cnstr.max_ineqRow_k, mpc.su_col, mpc.iS_0, mpc.iS_k, mpc.iS_ri_hat_0, mpc.iS_ri_hat_k);
end

function [ru_0, rse_k, ru_k, R_0, Q_k, R_k, Y_k] = reduced_kkt_control_rate_local(N, nu, R_0, Q_k, R_k, Y_k, ru_0, rse_k, ru_k, min_limit, max_limit, min_ineqRow_0, max_ineqRow_0, min_ineqRow_k, max_ineqRow_k, su_col, iS_0, iS_k, iS_ri_hat_0, iS_ri_hat_k)
if ~isempty(min_limit) && ~isempty(max_limit)

    % k = 0
    for i = 1:nu
        row_min_i = min_ineqRow_0(i);
        row_max_i = max_ineqRow_0(i);

        iS_min_i = iS_0(row_min_i);
        iS_max_i = iS_0(row_max_i);
        iS_ri_hat_min_i = iS_ri_hat_0(row_min_i);
        iS_ri_hat_max_i = iS_ri_hat_0(row_max_i);

        ru_0(i) = ru_0(i) + ...
            iS_ri_hat_max_i - iS_ri_hat_min_i;
        R_0(i,i) = R_0(i,i) + iS_min_i + iS_max_i;

    end

    % k = 1,...,N-1
    for k = 1:N-1
        for i = 1:nu
            su_i = su_col(i);
            row_min_i = min_ineqRow_k(i);
            row_max_i = max_ineqRow_k(i);

            iS_min_i = iS_k(row_min_i,k);
            iS_max_i = iS_k(row_max_i,k);
            iS_ri_hat_min_i = iS_ri_hat_k(row_min_i,k);
            iS_ri_hat_max_i = iS_ri_hat_k(row_max_i,k);

            iS_ri_hat_delta_i = iS_ri_hat_max_i - iS_ri_hat_min_i;
            iS_sum_i = iS_min_i + iS_max_i;

            rse_k(su_i,k) = rse_k(su_i,k) - ...
                iS_ri_hat_delta_i;
            ru_k(i,k) = ru_k(i,k) + iS_ri_hat_delta_i;

            Q_k(su_i,su_i,k) = Q_k(su_i,su_i,k) + iS_sum_i;
            R_k(i,i,k) = R_k(i,i,k) + iS_sum_i;
            Y_k(i,su_i,k) = Y_k(i,su_i,k) - iS_sum_i;

        end
    end
elseif ~isempty(min_limit)

    % k = 0
    for i = 1:nu
        row_i = min_ineqRow_0(i);

        iS_i = iS_0(row_i);
        iS_ri_hat_i = iS_ri_hat_0(row_i);

        % A_i = -I <-- u
        ru_0(i) = ru_0(i) - iS_ri_hat_i;
        R_0(i,i) =  R_0(i,i) + iS_i;

    end

    % k = 1,...,N-1
    for k = 1:N-1
        for i = 1:nu

            su_i = su_col(i);
            row_i = min_ineqRow_k(i);

            iS_i = iS_k(row_i,k);
            iS_ri_hat_i = iS_ri_hat_k(row_i,k);

            % A_i = [I -I] <-- [su u]^T
            rse_k(su_i,k) = rse_k(su_i,k) + iS_ri_hat_i;
            ru_k(i,k) = ru_k(i,k) - iS_ri_hat_i;

            Q_k(su_i,su_i,k) = Q_k(su_i,su_i,k) + iS_i;
            R_k(i,i,k) =  R_k(i,i,k) + iS_i;
            Y_k(i,su_i,k) = Y_k(i,su_i,k) - iS_i;

        end
    end
elseif ~isempty(max_limit)

    % k = 0
    for i = 1:nu
        row_i = max_ineqRow_0(i);

        iS_i = iS_0(row_i);
        iS_ri_hat_i = iS_ri_hat_0(row_i);

        % A_i = I <-- u
        ru_0(i) = ru_0(i) + iS_ri_hat_i;
        R_0(i,i) =  R_0(i,i) + iS_i;

    end

    % k = 1,...,N-1
    for k = 1:N-1
        for i = 1:nu

            su_i = su_col(i);
            row_i = max_ineqRow_k(i);

            iS_i = iS_k(row_i,k);
            iS_ri_hat_i = iS_ri_hat_k(row_i,k);

            % A_i = [-I I] <-- [su u]^T
            rse_k(su_i,k) = rse_k(su_i,k) - iS_ri_hat_i;
            ru_k(i,k) = ru_k(i,k) + iS_ri_hat_i;

            Q_k(su_i,su_i,k) = Q_k(su_i,su_i,k) + iS_i;
            R_k(i,i,k) =  R_k(i,i,k) + iS_i;
            Y_k(i,su_i,k) = Y_k(i,su_i,k) - iS_i;

        end
    end
end
end
