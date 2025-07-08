function feas = feas_check_s_prev(s_prev,mpc)

feas = 1;

% Check feasibility of initial state
if ~isempty(mpc.s_cnstr)
    mpc.s_cnstr = fi_box_fun(mpc.s_cnstr,s_prev,mpc.nx,mpc.nx,0);
    feas = fi_box_is_feasible(mpc.s_cnstr);
end


end