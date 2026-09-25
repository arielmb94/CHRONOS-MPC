function mpc = reduced_kkt_output_wrapper(mpc,y_cnstr)
[mpc.ru_0,mpc.rse_k,mpc.ru_k,mpc.rse_ter,mpc.R_0,mpc.Q_k,mpc.R_k,mpc.Y_k,mpc.Q_ter] = reduced_kkt_output_local(mpc.N, mpc.R_0, mpc.Q_k, mpc.R_k, mpc.Y_k, mpc.Q_ter, mpc.ru_0, mpc.rse_k, mpc.ru_k, mpc.rse_ter, mpc.s_col, y_cnstr.min_limit, y_cnstr.max_limit, y_cnstr.use_k0, y_cnstr.use_s, y_cnstr.use_u, y_cnstr.use_ter, y_cnstr.min_ineqRow_0, y_cnstr.max_ineqRow_0, y_cnstr.min_ineqRow_k, y_cnstr.max_ineqRow_k, y_cnstr.min_ineqRow_ter, y_cnstr.max_ineqRow_ter, mpc.iS_0, mpc.iS_k, mpc.iS_ter, mpc.iS_ri_hat_0, mpc.iS_ri_hat_k, mpc.iS_ri_hat_ter, mpc.D_0, mpc.C, mpc.D, mpc.C_ter);
end

function [ru_0, rse_k, ru_k, rse_ter, R_0, Q_k, R_k, Y_k, Q_ter] = reduced_kkt_output_local(N, R_0, Q_k, R_k, Y_k, Q_ter, ru_0, rse_k, ru_k, rse_ter, s_col, min_limit, max_limit, use_k0, use_s, use_u, use_ter, min_ineqRow_0, max_ineqRow_0, min_ineqRow_k, max_ineqRow_k, min_ineqRow_ter, max_ineqRow_ter, iS_0, iS_k, iS_ter, iS_ri_hat_0, iS_ri_hat_k, iS_ri_hat_ter, D_0, C, D, C_ter)
if ~isempty(min_limit) && ~isempty(max_limit)
    if ~isempty(use_k0)
        iS_ri_hat_y0 = iS_ri_hat_0(max_ineqRow_0) - iS_ri_hat_0(min_ineqRow_0);

        iS_y0 = iS_0(min_ineqRow_0) + iS_0(max_ineqRow_0);

        ru_0 = ru_0 + D_0' * iS_ri_hat_y0;
        R_0  = R_0  + D_0' * (iS_y0 .* D_0);
    end

    % k = 1,...,N-1
    for k = 1:N-1

        iS_ri_hat_yk = iS_ri_hat_k(max_ineqRow_k,k)-iS_ri_hat_k(min_ineqRow_k,k);

        iS_yk = iS_k(min_ineqRow_k,k) + iS_k(max_ineqRow_k,k);

        if ~isempty(use_s) && ~isempty(use_u)

            rse_k(s_col,k) = rse_k(s_col,k) + C(:,:,k)' * iS_ri_hat_yk;
            ru_k(:,k) = ru_k(:,k) + D(:,:,k)' * iS_ri_hat_yk;

            iS_C_k = iS_yk .* C(:,:,k);

            Q_k(s_col,s_col,k)  = Q_k(s_col,s_col,k)  + C(:,:,k)' * iS_C_k;
            R_k(:,:,k)  = R_k(:,:,k)  + D(:,:,k)' * (iS_yk .* D(:,:,k));
            Y_k(:,s_col,k)  = Y_k(:,s_col,k)  + D(:,:,k)' * iS_C_k;
        elseif ~isempty(use_s)

            rse_k(s_col,k) = rse_k(s_col,k) + C(:,:,k)' * iS_ri_hat_yk;

            Q_k(s_col,s_col,k)  = Q_k(s_col,s_col,k)  + C(:,:,k)' * (iS_yk .* C(:,:,k));
        elseif ~isempty(use_u)

            ru_k(:,k) = ru_k(:,k) + D(:,:,k)' * iS_ri_hat_yk;

            R_k(:,:,k)  = R_k(:,:,k)  + D(:,:,k)' * (iS_yk .* D(:,:,k));
        end
    end

    % k = N
    if ~isempty(use_ter)
        iS_ri_hat_yter = iS_ri_hat_ter(max_ineqRow_ter) - iS_ri_hat_ter(min_ineqRow_ter);

        iS_yter = iS_ter(min_ineqRow_ter) + iS_ter(max_ineqRow_ter);

        rse_ter(s_col) = rse_ter(s_col) + C_ter' * iS_ri_hat_yter;
        Q_ter(s_col,s_col)  = Q_ter(s_col,s_col) + C_ter' * (iS_yter .* C_ter);
    end

elseif ~isempty(min_limit)

    if ~isempty(use_k0)
        ru_0 = ru_0 - D_0' * iS_ri_hat_0(min_ineqRow_0);
        R_0  = R_0  + D_0' * (iS_0(min_ineqRow_0) .* D_0);
    end

    % k = 1,...,N-1
    for k = 1:N-1

        if ~isempty(use_s) && ~isempty(use_u)

            rse_k(s_col,k) = rse_k(s_col,k) - C(:,:,k)' * iS_ri_hat_k(min_ineqRow_k,k);
            ru_k(:,k) = ru_k(:,k) - D(:,:,k)' * iS_ri_hat_k(min_ineqRow_k,k);

            iS_C_k = iS_k(min_ineqRow_k,k) .* C(:,:,k);

            Q_k(s_col,s_col,k)  = Q_k(s_col,s_col,k) + C(:,:,k)' * iS_C_k;
            R_k(:,:,k)  = R_k(:,:,k) + D(:,:,k)' * (iS_k(min_ineqRow_k,k) .* D(:,:,k));
            Y_k(:,s_col,k)  = Y_k(:,s_col,k) + D(:,:,k)' * iS_C_k;
        elseif ~isempty(use_s)

            rse_k(s_col,k) = rse_k(s_col,k) - C(:,:,k)' * iS_ri_hat_k(min_ineqRow_k,k);

            Q_k(s_col,s_col,k)  = Q_k(s_col,s_col,k) + C(:,:,k)' * (iS_k(min_ineqRow_k,k) .* C(:,:,k));
        elseif ~isempty(use_u)

            ru_k(:,k) = ru_k(:,k) - D(:,:,k)' * iS_ri_hat_k(min_ineqRow_k,k);

            R_k(:,:,k) = R_k(:,:,k) + D(:,:,k)' * (iS_k(min_ineqRow_k,k) .* D(:,:,k));
        end
    end

    % k = N
    if ~isempty(use_ter)
        rse_ter(s_col) = rse_ter(s_col) - C_ter' * iS_ri_hat_ter(min_ineqRow_ter);
        Q_ter(s_col,s_col)  = Q_ter(s_col,s_col) + C_ter' * (iS_ter(min_ineqRow_ter) .* C_ter);
    end

elseif ~isempty(max_limit)

    if ~isempty(use_k0)
        ru_0 = ru_0 + D_0' * iS_ri_hat_0(max_ineqRow_0);
        R_0  = R_0  + D_0' * (iS_0(max_ineqRow_0) .* D_0);
    end

    % k = 1,...,N-1
    for k = 1:N-1

        if ~isempty(use_s) && ~isempty(use_u)

            rse_k(s_col,k) = rse_k(s_col,k) + C(:,:,k)' * iS_ri_hat_k(max_ineqRow_k,k);
            ru_k(:,k) = ru_k(:,k) + D(:,:,k)' * iS_ri_hat_k(max_ineqRow_k,k);

            iS_C_k = iS_k(max_ineqRow_k,k) .* C(:,:,k);

            Q_k(s_col,s_col,k)  = Q_k(s_col,s_col,k) + C(:,:,k)' * iS_C_k;
            R_k(:,:,k)  = R_k(:,:,k) + D(:,:,k)' * (iS_k(max_ineqRow_k,k) .* D(:,:,k));
            Y_k(:,s_col,k)  = Y_k(:,s_col,k) + D(:,:,k)' * iS_C_k;
        elseif ~isempty(use_s)

            rse_k(s_col,k) = rse_k(s_col,k) + C(:,:,k)' * iS_ri_hat_k(max_ineqRow_k,k);

            Q_k(s_col,s_col,k)  = Q_k(s_col,s_col,k) + C(:,:,k)' * (iS_k(max_ineqRow_k,k) .* C(:,:,k));
        elseif ~isempty(use_u)

            ru_k(:,k) = ru_k(:,k) + D(:,:,k)' * iS_ri_hat_k(max_ineqRow_k,k);

            R_k(:,:,k) = R_k(:,:,k) + D(:,:,k)' * (iS_k(max_ineqRow_k,k) .* D(:,:,k));
        end
    end

    % k = N
    if ~isempty(use_ter)
        rse_ter(s_col) = rse_ter(s_col) + C_ter' * iS_ri_hat_ter(max_ineqRow_ter);
        Q_ter(s_col,s_col)  = Q_ter(s_col,s_col) + C_ter' * (iS_ter(max_ineqRow_ter) .* C_ter);
    end
end
end
