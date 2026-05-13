function mpc = get_mpc_v(x,mpc)
    
    % extract updated slack values
    mpc.v(:) = x(mpc.Nx+mpc.Nu+1:mpc.Nx+mpc.Nu+mpc.Nv);

    if ~isempty(mpc.s_cnstr)
        if mpc.s_cnstr.min_limit
            mpc.s_cnstr = map_min_slack(mpc,mpc.s_cnstr);
        end
        if mpc.s_cnstr.max_limit
            mpc.s_cnstr = map_max_slack(mpc,mpc.s_cnstr);
        end
    end

    if ~isempty(mpc.y_cnstr)
        if mpc.y_cnstr.min_limit
            mpc.y_cnstr = map_min_slack(mpc,mpc.y_cnstr);
        end
        if mpc.y_cnstr.max_limit
            mpc.y_cnstr = map_max_slack(mpc,mpc.y_cnstr);
        end
    end

    if ~isempty(mpc.h_cnstr)
        if mpc.h_cnstr.min_limit
            mpc.h_cnstr = map_min_slack(mpc,mpc.h_cnstr);
        end
        if mpc.h_cnstr.max_limit
            mpc.h_cnstr = map_max_slack(mpc,mpc.h_cnstr);
        end
    end

    if mpc.ter_constraint
        mpc.v_ter = mpc.v(mpc.v_ter_global_index);
    end
    
end
