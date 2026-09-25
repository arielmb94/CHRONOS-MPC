function mpc = inequality_control_rate_wrapper(mpc,du_cnstr)
[mpc.ri_0,mpc.ri_k] = inequality_control_rate_local(mpc.ri_0,mpc.ri_k,...
    du_cnstr.min_limit,du_cnstr.max_limit,du_cnstr.min,du_cnstr.max,...
    du_cnstr.min_ineqRow_0,du_cnstr.min_ineqRow_k,...
    du_cnstr.max_ineqRow_0,du_cnstr.max_ineqRow_k,mpc.du,mpc.N);
end

function [ri_0,ri_k] = inequality_control_rate_local(ri_0,ri_k,min_limit,max_limit,min,max,...
    min_ineqRow_0,min_ineqRow_k,max_ineqRow_0,max_ineqRow_k,du,N)
if ~isempty(min_limit)
    ri_u = min - du;
    ri_0(min_ineqRow_0) = ri_0(min_ineqRow_0) + ri_u(:,1);
    ri_k(min_ineqRow_k,:) = ri_k(min_ineqRow_k,:) + ri_u(:,2:N);
end
if ~isempty(max_limit)
    ri_u = du-max;
    ri_0(max_ineqRow_0) = ri_0(max_ineqRow_0) + ri_u(:,1);
    ri_k(max_ineqRow_k,:) = ri_k(max_ineqRow_k,:) + ri_u(:,2:N);
end
end
