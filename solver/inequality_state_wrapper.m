function mpc = inequality_state_wrapper(mpc,s_cnstr)
[mpc.ri_k,mpc.ri_ter] = inequality_state_local(mpc.ri_k,mpc.ri_ter,...
    s_cnstr.min_limit,s_cnstr.max_limit,s_cnstr.min,s_cnstr.max,...
    s_cnstr.min_ineqRow_k,s_cnstr.min_ineqRow_ter,...
    s_cnstr.max_ineqRow_k,s_cnstr.max_ineqRow_ter,mpc.s,mpc.N);
end

function [ri_k,ri_ter] = inequality_state_local(ri_k,ri_ter,min_limit,max_limit,min,max,...
    min_ineqRow_k,min_ineqRow_ter,max_ineqRow_k,max_ineqRow_ter,s,N)
if ~isempty(min_limit)
    ri_s = min-s;
    ri_k(min_ineqRow_k,:) = ri_k(min_ineqRow_k,:) + ri_s(:,1:N-1);
    ri_ter(min_ineqRow_ter) = ri_ter(min_ineqRow_ter) + ri_s(:,N);
end
if ~isempty(max_limit)
    ri_s = s-max;
    ri_k(max_ineqRow_k,:) = ri_k(max_ineqRow_k,:) + ri_s(:,1:N-1);
    ri_ter(max_ineqRow_ter) = ri_ter(max_ineqRow_ter) + ri_s(:,N);
end
end
