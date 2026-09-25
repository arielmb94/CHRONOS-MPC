function mpc = inequality_custom_wrapper(mpc,h_cnstr)
[mpc.ri_0,mpc.ri_k,mpc.ri_ter] = inequality_custom_local(mpc.ri_0,mpc.ri_k,mpc.ri_ter,...
    h_cnstr.min_limit,h_cnstr.max_limit,...
    h_cnstr.min_0,h_cnstr.min,h_cnstr.min_ter,...
    h_cnstr.max_0,h_cnstr.max,h_cnstr.max_ter,...
    h_cnstr.min_ineqRow_0,h_cnstr.min_ineqRow_k,h_cnstr.min_ineqRow_ter,...
    h_cnstr.max_ineqRow_0,h_cnstr.max_ineqRow_k,h_cnstr.max_ineqRow_ter,...
    h_cnstr.use_k0,h_cnstr.use_ter,mpc.h_0,mpc.h,mpc.h_ter);
end

function [ri_0,ri_k,ri_ter] = inequality_custom_local(ri_0,ri_k,ri_ter,...
    min_limit,max_limit,min_0,min_k,min_ter,max_0,max_k,max_ter,...
    min_ineqRow_0,min_ineqRow_k,min_ineqRow_ter,...
    max_ineqRow_0,max_ineqRow_k,max_ineqRow_ter,use_k0,use_ter,h_0,h,h_ter)
if ~isempty(min_limit)
    if ~isempty(use_k0)
        ri_0(min_ineqRow_0) = ri_0(min_ineqRow_0) + min_0-h_0;
    end
    ri_k(min_ineqRow_k,:) = ri_k(min_ineqRow_k,:) + min_k-h;
    if ~isempty(use_ter)
        ri_ter(min_ineqRow_ter) = ri_ter(min_ineqRow_ter) + min_ter-h_ter;
    end
end
if ~isempty(max_limit)
    if ~isempty(use_k0)
        ri_0(max_ineqRow_0) = ri_0(max_ineqRow_0) + h_0-max_0;
    end
    ri_k(max_ineqRow_k,:) = ri_k(max_ineqRow_k,:) + h-max_k;
    if ~isempty(use_ter)
        ri_ter(max_ineqRow_ter) = ri_ter(max_ineqRow_ter) + h_ter-max_ter;
    end
end
end
