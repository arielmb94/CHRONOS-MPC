function mpc = update_mpc_state_cnstr(mpc,x_min,x_max)

if ~isempty(x_min)    
    mpc.x_min = x_min;
end

if ~isempty(x_max)    
    mpc.x_max = x_max;
end

end