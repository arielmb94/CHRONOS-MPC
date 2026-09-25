function mpc = reduced_kkt_custom_wrapper(mpc,h_cnstr)
[mpc.ru_0,mpc.rse_k,mpc.ru_k,mpc.rse_ter,mpc.R_0,mpc.Q_k,mpc.R_k,mpc.Y_k,mpc.Q_ter] = reduced_kkt_custom_local(mpc.N, mpc.R_0, mpc.Q_k, mpc.R_k, mpc.Y_k, mpc.Q_ter, mpc.ru_0, mpc.rse_k, mpc.ru_k, mpc.rse_ter, mpc.s_col, mpc.su_col, h_cnstr.min_limit, h_cnstr.max_limit, h_cnstr.use_k0, h_cnstr.use_s, h_cnstr.use_su, h_cnstr.use_u, h_cnstr.use_ter, h_cnstr.min_ineqRow_0, h_cnstr.max_ineqRow_0, h_cnstr.min_ineqRow_k, h_cnstr.max_ineqRow_k, h_cnstr.min_ineqRow_ter, h_cnstr.max_ineqRow_ter, mpc.iS_0, mpc.iS_k, mpc.iS_ter, mpc.iS_ri_hat_0, mpc.iS_ri_hat_k, mpc.iS_ri_hat_ter, mpc.Dh_0, mpc.Ch, mpc.Dsuh, mpc.Dh, mpc.Ch_ter);
end

function [ru_0, rse_k, ru_k, rse_ter, R_0, Q_k, R_k, Y_k, Q_ter] = reduced_kkt_custom_local(N, R_0, Q_k, R_k, Y_k, Q_ter, ru_0, rse_k, ru_k, rse_ter, s_col, su_col, min_limit, max_limit, use_k0, use_s, use_su, use_u, use_ter, min_ineqRow_0, max_ineqRow_0, min_ineqRow_k, max_ineqRow_k, min_ineqRow_ter, max_ineqRow_ter, iS_0, iS_k, iS_ter, iS_ri_hat_0, iS_ri_hat_k, iS_ri_hat_ter, Dh_0, Ch, Dsuh, Dh, Ch_ter)
if ~isempty(min_limit) && ~isempty(max_limit)
    if ~isempty(use_k0)
        iS_ri_hat_h0 = iS_ri_hat_0(max_ineqRow_0) - iS_ri_hat_0(min_ineqRow_0);
        iS_h0 = iS_0(min_ineqRow_0) + iS_0(max_ineqRow_0);

        ru_0 = ru_0 + Dh_0' * iS_ri_hat_h0;
        R_0 = R_0 + Dh_0' * (iS_h0 .* Dh_0);
    end

    % k = 1,...,N-1
    for k = 1:N-1

        iS_ri_hat_hk = iS_ri_hat_k(max_ineqRow_k,k) - iS_ri_hat_k(min_ineqRow_k,k);
        iS_hk = iS_k(min_ineqRow_k,k) + iS_k(max_ineqRow_k,k);

        if ~isempty(use_s) && ~isempty(use_su) && ~isempty(use_u)

            rse_k(s_col,k) = rse_k(s_col,k) + Ch(:,:,k)' * iS_ri_hat_hk;
            rse_k(su_col,k) = rse_k(su_col,k) + Dsuh(:,:,k)' * iS_ri_hat_hk;
            ru_k(:,k) = ru_k(:,k) + Dh(:,:,k)' * iS_ri_hat_hk;

            iS_Ch_k = iS_hk .* Ch(:,:,k);
            iS_Dsuh_k = iS_hk .* Dsuh(:,:,k);

            Q_k(s_col,s_col,k) = Q_k(s_col,s_col,k) + Ch(:,:,k)' * iS_Ch_k;
            Q_k(s_col,su_col,k) = Q_k(s_col,su_col,k) + Ch(:,:,k)' * iS_Dsuh_k;
            Q_k(su_col,s_col,k) = Q_k(su_col,s_col,k) + Dsuh(:,:,k)' * iS_Ch_k;
            Q_k(su_col,su_col,k) = Q_k(su_col,su_col,k) + Dsuh(:,:,k)' * iS_Dsuh_k;
            R_k(:,:,k) = R_k(:,:,k) + Dh(:,:,k)' * (iS_hk .* Dh(:,:,k));
            Y_k(:,s_col,k) = Y_k(:,s_col,k) + Dh(:,:,k)' * iS_Ch_k;
            Y_k(:,su_col,k) = Y_k(:,su_col,k) + Dh(:,:,k)' * iS_Dsuh_k;
        elseif ~isempty(use_s) && ~isempty(use_su)

            rse_k(s_col,k) = rse_k(s_col,k) + Ch(:,:,k)' * iS_ri_hat_hk;
            rse_k(su_col,k) = rse_k(su_col,k) + Dsuh(:,:,k)' * iS_ri_hat_hk;

            iS_Ch_k = iS_hk .* Ch(:,:,k);
            iS_Dsuh_k = iS_hk .* Dsuh(:,:,k);

            Q_k(s_col,s_col,k) = Q_k(s_col,s_col,k) + Ch(:,:,k)' * iS_Ch_k;
            Q_k(s_col,su_col,k) = Q_k(s_col,su_col,k) + Ch(:,:,k)' * iS_Dsuh_k;
            Q_k(su_col,s_col,k) = Q_k(su_col,s_col,k) + Dsuh(:,:,k)' * iS_Ch_k;
            Q_k(su_col,su_col,k) = Q_k(su_col,su_col,k) + Dsuh(:,:,k)' * iS_Dsuh_k;
        elseif ~isempty(use_s) && ~isempty(use_u)

            rse_k(s_col,k) = rse_k(s_col,k) + Ch(:,:,k)' * iS_ri_hat_hk;
            ru_k(:,k) = ru_k(:,k) + Dh(:,:,k)' * iS_ri_hat_hk;

            iS_Ch_k = iS_hk .* Ch(:,:,k);

            Q_k(s_col,s_col,k) = Q_k(s_col,s_col,k) + Ch(:,:,k)' * iS_Ch_k;
            R_k(:,:,k) = R_k(:,:,k) + Dh(:,:,k)' * (iS_hk .* Dh(:,:,k));
            Y_k(:,s_col,k) = Y_k(:,s_col,k) + Dh(:,:,k)' * iS_Ch_k;
        elseif ~isempty(use_su) && ~isempty(use_u)

            rse_k(su_col,k) = rse_k(su_col,k) + Dsuh(:,:,k)' * iS_ri_hat_hk;
            ru_k(:,k) = ru_k(:,k) + Dh(:,:,k)' * iS_ri_hat_hk;

            iS_Dsuh_k = iS_hk .* Dsuh(:,:,k);

            Q_k(su_col,su_col,k) = Q_k(su_col,su_col,k) + Dsuh(:,:,k)' * iS_Dsuh_k;
            R_k(:,:,k) = R_k(:,:,k) + Dh(:,:,k)' * (iS_hk .* Dh(:,:,k));
            Y_k(:,su_col,k) = Y_k(:,su_col,k) + Dh(:,:,k)' * iS_Dsuh_k;
        elseif ~isempty(use_s)

            rse_k(s_col,k) = rse_k(s_col,k) + Ch(:,:,k)' * iS_ri_hat_hk;
            Q_k(s_col,s_col,k) = Q_k(s_col,s_col,k) + Ch(:,:,k)' * (iS_hk .* Ch(:,:,k));
        elseif ~isempty(use_su)

            rse_k(su_col,k) = rse_k(su_col,k) + Dsuh(:,:,k)' * iS_ri_hat_hk;
            Q_k(su_col,su_col,k) = Q_k(su_col,su_col,k) + Dsuh(:,:,k)' * (iS_hk .* Dsuh(:,:,k));
        elseif ~isempty(use_u)

            ru_k(:,k) = ru_k(:,k) + Dh(:,:,k)' * iS_ri_hat_hk;
            R_k(:,:,k) = R_k(:,:,k) + Dh(:,:,k)' * (iS_hk .* Dh(:,:,k));
        end
    end

    % k = N
    if ~isempty(use_ter)
        iS_ri_hat_hter = iS_ri_hat_ter(max_ineqRow_ter) - iS_ri_hat_ter(min_ineqRow_ter);
        iS_hter = iS_ter(min_ineqRow_ter) + iS_ter(max_ineqRow_ter);

        rse_ter(s_col) = rse_ter(s_col) + Ch_ter' * iS_ri_hat_hter;
        Q_ter(s_col,s_col) = Q_ter(s_col,s_col) + Ch_ter' * (iS_hter .* Ch_ter);
    end

elseif ~isempty(min_limit)

    if ~isempty(use_k0)
        ru_0 = ru_0 - Dh_0' * iS_ri_hat_0(min_ineqRow_0);
        R_0 = R_0 + Dh_0' * (iS_0(min_ineqRow_0) .* Dh_0);
    end

    % k = 1,...,N-1
    for k = 1:N-1

        if ~isempty(use_s) && ~isempty(use_su) && ~isempty(use_u)

            rse_k(s_col,k) = rse_k(s_col,k) - Ch(:,:,k)' * iS_ri_hat_k(min_ineqRow_k,k);
            rse_k(su_col,k) = rse_k(su_col,k) - Dsuh(:,:,k)' * iS_ri_hat_k(min_ineqRow_k,k);
            ru_k(:,k) = ru_k(:,k) - Dh(:,:,k)' * iS_ri_hat_k(min_ineqRow_k,k);

            iS_Ch_k = iS_k(min_ineqRow_k,k) .* Ch(:,:,k);
            iS_Dsuh_k = iS_k(min_ineqRow_k,k) .* Dsuh(:,:,k);

            Q_k(s_col,s_col,k) = Q_k(s_col,s_col,k) + Ch(:,:,k)' * iS_Ch_k;
            Q_k(s_col,su_col,k) = Q_k(s_col,su_col,k) + Ch(:,:,k)' * iS_Dsuh_k;
            Q_k(su_col,s_col,k) = Q_k(su_col,s_col,k) + Dsuh(:,:,k)' * iS_Ch_k;
            Q_k(su_col,su_col,k) = Q_k(su_col,su_col,k) + Dsuh(:,:,k)' * iS_Dsuh_k;
            R_k(:,:,k) = R_k(:,:,k) + Dh(:,:,k)' * (iS_k(min_ineqRow_k,k) .* Dh(:,:,k));
            Y_k(:,s_col,k) = Y_k(:,s_col,k) + Dh(:,:,k)' * iS_Ch_k;
            Y_k(:,su_col,k) = Y_k(:,su_col,k) + Dh(:,:,k)' * iS_Dsuh_k;
        elseif ~isempty(use_s) && ~isempty(use_su)

            rse_k(s_col,k) = rse_k(s_col,k) - Ch(:,:,k)' * iS_ri_hat_k(min_ineqRow_k,k);
            rse_k(su_col,k) = rse_k(su_col,k) - Dsuh(:,:,k)' * iS_ri_hat_k(min_ineqRow_k,k);

            iS_Ch_k = iS_k(min_ineqRow_k,k) .* Ch(:,:,k);
            iS_Dsuh_k = iS_k(min_ineqRow_k,k) .* Dsuh(:,:,k);

            Q_k(s_col,s_col,k) = Q_k(s_col,s_col,k) + Ch(:,:,k)' * iS_Ch_k;
            Q_k(s_col,su_col,k) = Q_k(s_col,su_col,k) + Ch(:,:,k)' * iS_Dsuh_k;
            Q_k(su_col,s_col,k) = Q_k(su_col,s_col,k) + Dsuh(:,:,k)' * iS_Ch_k;
            Q_k(su_col,su_col,k) = Q_k(su_col,su_col,k) + Dsuh(:,:,k)' * iS_Dsuh_k;
        elseif ~isempty(use_s) && ~isempty(use_u)

            rse_k(s_col,k) = rse_k(s_col,k) - Ch(:,:,k)' * iS_ri_hat_k(min_ineqRow_k,k);
            ru_k(:,k) = ru_k(:,k) - Dh(:,:,k)' * iS_ri_hat_k(min_ineqRow_k,k);

            iS_Ch_k = iS_k(min_ineqRow_k,k) .* Ch(:,:,k);

            Q_k(s_col,s_col,k) = Q_k(s_col,s_col,k) + Ch(:,:,k)' * iS_Ch_k;
            R_k(:,:,k) = R_k(:,:,k) + Dh(:,:,k)' * (iS_k(min_ineqRow_k,k) .* Dh(:,:,k));
            Y_k(:,s_col,k) = Y_k(:,s_col,k) + Dh(:,:,k)' * iS_Ch_k;
        elseif ~isempty(use_su) && ~isempty(use_u)

            rse_k(su_col,k) = rse_k(su_col,k) - Dsuh(:,:,k)' * iS_ri_hat_k(min_ineqRow_k,k);
            ru_k(:,k) = ru_k(:,k) - Dh(:,:,k)' * iS_ri_hat_k(min_ineqRow_k,k);

            iS_Dsuh_k = iS_k(min_ineqRow_k,k) .* Dsuh(:,:,k);

            Q_k(su_col,su_col,k) = Q_k(su_col,su_col,k) + Dsuh(:,:,k)' * iS_Dsuh_k;
            R_k(:,:,k) = R_k(:,:,k) + Dh(:,:,k)' * (iS_k(min_ineqRow_k,k) .* Dh(:,:,k));
            Y_k(:,su_col,k) = Y_k(:,su_col,k) + Dh(:,:,k)' * iS_Dsuh_k;
        elseif ~isempty(use_s)

            rse_k(s_col,k) = rse_k(s_col,k) - Ch(:,:,k)' * iS_ri_hat_k(min_ineqRow_k,k);
            Q_k(s_col,s_col,k) = Q_k(s_col,s_col,k) + Ch(:,:,k)' * (iS_k(min_ineqRow_k,k) .* Ch(:,:,k));
        elseif ~isempty(use_su)

            rse_k(su_col,k) = rse_k(su_col,k) - Dsuh(:,:,k)' * iS_ri_hat_k(min_ineqRow_k,k);
            Q_k(su_col,su_col,k) = Q_k(su_col,su_col,k) + Dsuh(:,:,k)' * (iS_k(min_ineqRow_k,k) .* Dsuh(:,:,k));
        elseif ~isempty(use_u)

            ru_k(:,k) = ru_k(:,k) - Dh(:,:,k)' * iS_ri_hat_k(min_ineqRow_k,k);
            R_k(:,:,k) = R_k(:,:,k) + Dh(:,:,k)' * (iS_k(min_ineqRow_k,k) .* Dh(:,:,k));
        end
    end

    % k = N
    if ~isempty(use_ter)
        rse_ter(s_col) = rse_ter(s_col) - Ch_ter' * iS_ri_hat_ter(min_ineqRow_ter);
        Q_ter(s_col,s_col) = Q_ter(s_col,s_col) + Ch_ter' * (iS_ter(min_ineqRow_ter) .* Ch_ter);
    end

elseif ~isempty(max_limit)

    if ~isempty(use_k0)
        ru_0 = ru_0 + Dh_0' * iS_ri_hat_0(max_ineqRow_0);
        R_0 = R_0 + Dh_0' * (iS_0(max_ineqRow_0) .* Dh_0);
    end

    % k = 1,...,N-1
    for k = 1:N-1

        if ~isempty(use_s) && ~isempty(use_su) && ~isempty(use_u)

            rse_k(s_col,k) = rse_k(s_col,k) + Ch(:,:,k)' * iS_ri_hat_k(max_ineqRow_k,k);
            rse_k(su_col,k) = rse_k(su_col,k) + Dsuh(:,:,k)' * iS_ri_hat_k(max_ineqRow_k,k);
            ru_k(:,k) = ru_k(:,k) + Dh(:,:,k)' * iS_ri_hat_k(max_ineqRow_k,k);

            iS_Ch_k = iS_k(max_ineqRow_k,k) .* Ch(:,:,k);
            iS_Dsuh_k = iS_k(max_ineqRow_k,k) .* Dsuh(:,:,k);

            Q_k(s_col,s_col,k) = Q_k(s_col,s_col,k) + Ch(:,:,k)' * iS_Ch_k;
            Q_k(s_col,su_col,k) = Q_k(s_col,su_col,k) + Ch(:,:,k)' * iS_Dsuh_k;
            Q_k(su_col,s_col,k) = Q_k(su_col,s_col,k) + Dsuh(:,:,k)' * iS_Ch_k;
            Q_k(su_col,su_col,k) = Q_k(su_col,su_col,k) + Dsuh(:,:,k)' * iS_Dsuh_k;
            R_k(:,:,k) = R_k(:,:,k) + Dh(:,:,k)' * (iS_k(max_ineqRow_k,k) .* Dh(:,:,k));
            Y_k(:,s_col,k) = Y_k(:,s_col,k) + Dh(:,:,k)' * iS_Ch_k;
            Y_k(:,su_col,k) = Y_k(:,su_col,k) + Dh(:,:,k)' * iS_Dsuh_k;
        elseif ~isempty(use_s) && ~isempty(use_su)

            rse_k(s_col,k) = rse_k(s_col,k) + Ch(:,:,k)' * iS_ri_hat_k(max_ineqRow_k,k);
            rse_k(su_col,k) = rse_k(su_col,k) + Dsuh(:,:,k)' * iS_ri_hat_k(max_ineqRow_k,k);

            iS_Ch_k = iS_k(max_ineqRow_k,k) .* Ch(:,:,k);
            iS_Dsuh_k = iS_k(max_ineqRow_k,k) .* Dsuh(:,:,k);

            Q_k(s_col,s_col,k) = Q_k(s_col,s_col,k) + Ch(:,:,k)' * iS_Ch_k;
            Q_k(s_col,su_col,k) = Q_k(s_col,su_col,k) + Ch(:,:,k)' * iS_Dsuh_k;
            Q_k(su_col,s_col,k) = Q_k(su_col,s_col,k) + Dsuh(:,:,k)' * iS_Ch_k;
            Q_k(su_col,su_col,k) = Q_k(su_col,su_col,k) + Dsuh(:,:,k)' * iS_Dsuh_k;
        elseif ~isempty(use_s) && ~isempty(use_u)

            rse_k(s_col,k) = rse_k(s_col,k) + Ch(:,:,k)' * iS_ri_hat_k(max_ineqRow_k,k);
            ru_k(:,k) = ru_k(:,k) + Dh(:,:,k)' * iS_ri_hat_k(max_ineqRow_k,k);

            iS_Ch_k = iS_k(max_ineqRow_k,k) .* Ch(:,:,k);

            Q_k(s_col,s_col,k) = Q_k(s_col,s_col,k) + Ch(:,:,k)' * iS_Ch_k;
            R_k(:,:,k) = R_k(:,:,k) + Dh(:,:,k)' * (iS_k(max_ineqRow_k,k) .* Dh(:,:,k));
            Y_k(:,s_col,k) = Y_k(:,s_col,k) + Dh(:,:,k)' * iS_Ch_k;
        elseif ~isempty(use_su) && ~isempty(use_u)

            rse_k(su_col,k) = rse_k(su_col,k) + Dsuh(:,:,k)' * iS_ri_hat_k(max_ineqRow_k,k);
            ru_k(:,k) = ru_k(:,k) + Dh(:,:,k)' * iS_ri_hat_k(max_ineqRow_k,k);

            iS_Dsuh_k = iS_k(max_ineqRow_k,k) .* Dsuh(:,:,k);

            Q_k(su_col,su_col,k) = Q_k(su_col,su_col,k) + Dsuh(:,:,k)' * iS_Dsuh_k;
            R_k(:,:,k) = R_k(:,:,k) + Dh(:,:,k)' * (iS_k(max_ineqRow_k,k) .* Dh(:,:,k));
            Y_k(:,su_col,k) = Y_k(:,su_col,k) + Dh(:,:,k)' * iS_Dsuh_k;
        elseif ~isempty(use_s)

            rse_k(s_col,k) = rse_k(s_col,k) + Ch(:,:,k)' * iS_ri_hat_k(max_ineqRow_k,k);
            Q_k(s_col,s_col,k) = Q_k(s_col,s_col,k) + Ch(:,:,k)' * (iS_k(max_ineqRow_k,k) .* Ch(:,:,k));
        elseif ~isempty(use_su)

            rse_k(su_col,k) = rse_k(su_col,k) + Dsuh(:,:,k)' * iS_ri_hat_k(max_ineqRow_k,k);
            Q_k(su_col,su_col,k) = Q_k(su_col,su_col,k) + Dsuh(:,:,k)' * (iS_k(max_ineqRow_k,k) .* Dsuh(:,:,k));
        elseif ~isempty(use_u)

            ru_k(:,k) = ru_k(:,k) + Dh(:,:,k)' * iS_ri_hat_k(max_ineqRow_k,k);
            R_k(:,:,k) = R_k(:,:,k) + Dh(:,:,k)' * (iS_k(max_ineqRow_k,k) .* Dh(:,:,k));
        end
    end

    % k = N
    if ~isempty(use_ter)
        rse_ter(s_col) = rse_ter(s_col) + Ch_ter' * iS_ri_hat_ter(max_ineqRow_ter);
        Q_ter(s_col,s_col) = Q_ter(s_col,s_col) + Ch_ter' * (iS_ter(max_ineqRow_ter) .* Ch_ter);
    end
end
end
