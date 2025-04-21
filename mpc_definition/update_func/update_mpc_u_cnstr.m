function mpc = update_mpc_u_cnstr(mpc,u_min,u_max)

if ~isempty(u_min)    
    mpc.u_min = u_min;
end

if ~isempty(u_max)    
    mpc.u_max = u_max;
end

end