function mpc = reduced_kkt_state_wrapper(mpc,x_cnstr)
[mpc.rxe_k,mpc.rxe_ter,mpc.Q_k,mpc.Q_ter] = reduced_kkt_state_local(mpc.N, mpc.nx, mpc.Q_k, mpc.Q_ter, mpc.rxe_k, mpc.rxe_ter, x_cnstr.min_limit, x_cnstr.max_limit, x_cnstr.min_ineqRow_k, x_cnstr.max_ineqRow_k, x_cnstr.min_ineqRow_ter, x_cnstr.max_ineqRow_ter, mpc.x_col, mpc.iS_k, mpc.iS_ter, mpc.iS_ri_hat_k, mpc.iS_ri_hat_ter);
end

function [rxe_k, rxe_ter, Q_k, Q_ter] = reduced_kkt_state_local(N, nx, Q_k, Q_ter, rxe_k, rxe_ter, min_limit, max_limit, min_ineqRow_k, max_ineqRow_k, min_ineqRow_ter, max_ineqRow_ter, x_col, iS_k, iS_ter, iS_ri_hat_k, iS_ri_hat_ter)
if ~isempty(min_limit) && ~isempty(max_limit)

    % k = 1,...,N-1
    for k = 1:N-1
        for i = 1:nx
            x_i = x_col(i);
            row_min_i = min_ineqRow_k(i);
            row_max_i = max_ineqRow_k(i);

            iS_min_i = iS_k(row_min_i,k);
            iS_max_i = iS_k(row_max_i,k);
            iS_ri_hat_min_i = iS_ri_hat_k(row_min_i,k);
            iS_ri_hat_max_i = iS_ri_hat_k(row_max_i,k);

            rxe_k(x_i,k) = rxe_k(x_i,k) + ...
                iS_ri_hat_max_i - iS_ri_hat_min_i;
            Q_k(x_i,x_i,k) = Q_k(x_i,x_i,k) + ...
                iS_min_i + iS_max_i;

        end
    end

    % k = N
    for i = 1:nx
        row_min_i = min_ineqRow_ter(i);
        row_max_i = max_ineqRow_ter(i);

        iS_min_i = iS_ter(row_min_i);
        iS_max_i = iS_ter(row_max_i);
        iS_ri_hat_min_i = iS_ri_hat_ter(row_min_i);
        iS_ri_hat_max_i = iS_ri_hat_ter(row_max_i);

        rxe_ter(i) = rxe_ter(i) + ...
            iS_ri_hat_max_i - iS_ri_hat_min_i;
        Q_ter(i,i) = Q_ter(i,i) + iS_min_i + iS_max_i;

    end
elseif ~isempty(min_limit)

    % k = 1,...,N-1
    for k = 1:N-1
        for i = 1:nx

            x_i = x_col(i);
            row_i = min_ineqRow_k(i);

            iS_i = iS_k(row_i,k);
            iS_ri_hat_i = iS_ri_hat_k(row_i,k);

            % A_i = -I <-- x
            rxe_k(x_i,k) = rxe_k(x_i,k) - iS_ri_hat_i;

            Q_k(x_i,x_i,k) = Q_k(x_i,x_i,k) + iS_i;

        end
    end

    % k = N
    for i = 1:nx
        row_i = min_ineqRow_ter(i);

        iS_i = iS_ter(row_i);
        iS_ri_hat_i = iS_ri_hat_ter(row_i);

        % A_i = -I <-- x
        rxe_ter(i) = rxe_ter(i) - iS_ri_hat_i;
        Q_ter(i,i) =  Q_ter(i,i) + iS_i;

    end
elseif ~isempty(max_limit)

    % k = 1,...,N-1
    for k = 1:N-1
        for i = 1:nx

            x_i = x_col(i);
            row_i = max_ineqRow_k(i);

            iS_i = iS_k(row_i,k);
            iS_ri_hat_i = iS_ri_hat_k(row_i,k);

            % A_i = I <-- x
            rxe_k(x_i,k) = rxe_k(x_i,k) + iS_ri_hat_i;

            Q_k(x_i,x_i,k) = Q_k(x_i,x_i,k) + iS_i;

        end
    end

    % k = N
    for i = 1:nx
        row_i = max_ineqRow_ter(i);

        iS_i = iS_ter(row_i);
        iS_ri_hat_i = iS_ri_hat_ter(row_i);

        % A_i = I <-- x
        rxe_ter(i) = rxe_ter(i) + iS_ri_hat_i;
        Q_ter(i,i) =  Q_ter(i,i) + iS_i;

    end
end
end
