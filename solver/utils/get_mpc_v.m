function mpc = get_mpc_v(x,mpc)
    
    % extract updated slack values
    mpc.v = x(mpc.Nx+mpc.Nu+1:mpc.Nx+mpc.Nu+mpc.Nv);

    if ~isempty(mpc.s_cnstr)
        mpc.s_cnstr = map_slack(mpc,mpc.s_cnstr);
    end
    if ~isempty(mpc.s_ter_cnstr)
        mpc.s_ter_cnstr = map_slack(mpc,mpc.s_ter_cnstr);
    end
    if ~isempty(mpc.u_cnstr)
        mpc.u_cnstr = map_slack(mpc,mpc.u_cnstr);
    end
    if ~isempty(mpc.du_cnstr)
        mpc.du_cnstr = map_slack(mpc,mpc.du_cnstr);
    end
    if ~isempty(mpc.y_cnstr)
        mpc.y_cnstr = map_slack(mpc,mpc.y_cnstr);
    end
    if ~isempty(mpc.h_cnstr)
        mpc.h_cnstr = map_slack(mpc,mpc.h_cnstr);
    end
    if mpc.ter_constraint
        mpc.v_ter = mpc.ter_cnstr_map*mpc.v;
    end
    
end

% map slack value to constraint for box constraints
function cnstr = map_slack(mpc,cnstr)
    if cnstr.min_slack_nv
        cnstr.min_v = cnstr.min_slack_map'*mpc.v;
    end
    if cnstr.max_slack_nv
        cnstr.max_v = cnstr.max_slack_map'*mpc.v;
    end
end