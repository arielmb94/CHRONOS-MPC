function mpc = set_soft_cost_qv(mpc)

if ~isempty(mpc.has_s_cnstr)
    mpc.s_cnstr = init_gradSlack(mpc,mpc.s_cnstr);
end

if ~isempty(mpc.has_y_cnstr)
    mpc.y_cnstr = init_gradSlack(mpc,mpc.y_cnstr);
end

if ~isempty(mpc.has_h_cnstr)
    mpc.h_cnstr = init_gradSlack(mpc,mpc.h_cnstr);
end

end


function cnstr = init_gradSlack(mpc,cnstr)

if isempty(mpc.tracking_cost)
    qv = mpc.qv;
elseif mpc.N > 1
    max_eig_Qe = -inf;
    for k = 1:mpc.N-1
        max_eig_Qe = max(max_eig_Qe, max(real(eig(mpc.Qe(:,:,k)))));
    end
    qv = min(mpc.qv, max_eig_Qe);
elseif ~isempty(mpc.y_use_ter)
    max_eig_Qe_ter = max(real(eig(mpc.Qe_ter)));
    qv = min(mpc.qv, max_eig_Qe_ter);
elseif ~isempty(mpc.y_use_k0)
    max_eig_Qe_0 = max(real(eig(mpc.Qe_0)));
    qv = min(mpc.qv, max_eig_Qe_0);
else
    qv = mpc.qv;
end

if ~isempty(cnstr.min_limit)
    if ~isempty(cnstr.use_k0)
        if ~any(cnstr.qv_min_0(:))
            cnstr.qv_min_0(:) = qv;
        end
    end
    if ~any(cnstr.qv_min(:))
        cnstr.qv_min(:,:) = qv;
    end
    if ~isempty(cnstr.use_ter)
        if ~any(cnstr.qv_min_ter(:))
            cnstr.qv_min_ter(:) = qv;
        end
    end
end

if ~isempty(cnstr.max_limit)
    if ~isempty(cnstr.use_k0)
        if ~any(cnstr.qv_max_0(:))
            cnstr.qv_max_0(:) = qv;
        end
    end
    if ~any(cnstr.qv_max(:))
        cnstr.qv_max(:,:) = qv;
    end
    if ~isempty(cnstr.use_ter)
        if ~any(cnstr.qv_max_ter(:))
            cnstr.qv_max_ter(:) = qv;
        end
    end
end

end
