function feas = feas_check_s_prev(s_prev,mpc)

feas = 1;

% Check feasibility of initial state
if ~isempty(mpc.x_min) || ~isempty(mpc.x_max)
    [fi_s_prev_min_x0,fi_s_prev_max_x0] = fi_box_fun(s_prev,...
                                        mpc.x_min,mpc.x_max,mpc.nx,mpc.nx);
    if any(fi_s_prev_min_x0>0) || any(fi_s_prev_max_x0>0)
        feas = 0;
    end
end


end