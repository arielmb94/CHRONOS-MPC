function mpc = update_mpc_ter_state_cnstr(mpc,x_ter_min,x_ter_max)

if ~isempty(x_ter_min)    
    mpc.x_ter_min = x_ter_min;
end

if ~isempty(x_ter_max)    
    mpc.x_ter_max = x_ter_max;
end

end