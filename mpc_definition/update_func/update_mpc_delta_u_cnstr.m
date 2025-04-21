function mpc = update_mpc_delta_u_cnstr(mpc,du_min,du_max)

if ~isempty(du_min)    
    mpc.du_min = du_min;
end

if ~isempty(du_max)    
    mpc.du_max = du_max;
end

end