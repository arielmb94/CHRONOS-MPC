function mpc = update_mpc_output_cnstr(mpc,y_min,y_max)

if ~isempty(y_min)    
    mpc.y_min = y_min;
end

if ~isempty(y_max)    
    mpc.y_max = y_max;
end

end