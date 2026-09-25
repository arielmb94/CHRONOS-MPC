function mpc = inequality_control_wrapper(mpc,u_cnstr)
[mpc.ri_0,mpc.ri_k] = inequality_control_local(mpc.ri_0,mpc.ri_k,...
    u_cnstr.min_limit,u_cnstr.max_limit,u_cnstr.min,u_cnstr.max,...
    u_cnstr.min_ineqRow_0,u_cnstr.min_ineqRow_k,...
    u_cnstr.max_ineqRow_0,u_cnstr.max_ineqRow_k,mpc.u,mpc.N);
end

function [ri_0,ri_k] = inequality_control_local(ri_0,ri_k,min_limit,max_limit,min,max,...
    min_ineqRow_0,min_ineqRow_k,max_ineqRow_0,max_ineqRow_k,u,N)
if ~isempty(min_limit)
    ri_u = min - u;
    ri_0(min_ineqRow_0) = ri_0(min_ineqRow_0) + ri_u(:,1);
    ri_k(min_ineqRow_k,:) = ri_k(min_ineqRow_k,:) + ri_u(:,2:N);
end
if ~isempty(max_limit)
    ri_u = u-max;
    ri_0(max_ineqRow_0) = ri_0(max_ineqRow_0) + ri_u(:,1);
    ri_k(max_ineqRow_k,:) = ri_k(max_ineqRow_k,:) + ri_u(:,2:N);
end
end
