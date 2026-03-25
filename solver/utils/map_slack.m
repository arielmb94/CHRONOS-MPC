% map slack value to constraint for box constraints
function cnstr = map_slack(mpc,cnstr)
    if cnstr.min_slack_nv
        vi = mpc.v(cnstr.min_v_global_index);
        cnstr.min_v = cnstr.min_slack_local_map'*vi;
    end
    if cnstr.max_slack_nv
        vi = mpc.v(cnstr.max_v_global_index);
        cnstr.max_v = cnstr.max_slack_local_map'*vi;
    end
end