% map slack value to constraint for box constraints
function cnstr = map_max_slack(mpc,cnstr)

cnstr.max_v = mpc.v(cnstr.max_v_global_index);    

end