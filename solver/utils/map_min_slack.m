% map slack value to constraint for box constraints
function cnstr = map_min_slack(mpc,cnstr)

cnstr.min_v = mpc.v(cnstr.min_v_global_index);    

end