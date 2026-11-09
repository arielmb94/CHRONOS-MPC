function mpc = inequality_state_wrapper(mpc,x_cnstr)
[mpc.ri_k,mpc.ri_ter] = inequality_state_local(mpc.ri_k,mpc.ri_ter,...
    x_cnstr.min_limit,x_cnstr.max_limit,x_cnstr.min,x_cnstr.max,...
    x_cnstr.min_ineqRow_k,x_cnstr.min_ineqRow_ter,...
    x_cnstr.max_ineqRow_k,x_cnstr.max_ineqRow_ter,mpc.x,mpc.N);
end

function [ri_k,ri_ter] = inequality_state_local(ri_k,ri_ter,min_limit,max_limit,min,max,...
    min_ineqRow_k,min_ineqRow_ter,max_ineqRow_k,max_ineqRow_ter,x,N)
if ~isempty(min_limit)
    ri_x = min-x;
    ri_k(min_ineqRow_k,:) = ri_k(min_ineqRow_k,:) + ri_x(:,1:N-1);
    ri_ter(min_ineqRow_ter) = ri_ter(min_ineqRow_ter) + ri_x(:,N);
end
if ~isempty(max_limit)
    ri_x = x-max;
    ri_k(max_ineqRow_k,:) = ri_k(max_ineqRow_k,:) + ri_x(:,1:N-1);
    ri_ter(max_ineqRow_ter) = ri_ter(max_ineqRow_ter) + ri_x(:,N);
end
end
