function mpc = inequality_output_wrapper(mpc,y_cnstr)
[mpc.ri_0,mpc.ri_k,mpc.ri_ter] = inequality_output_local(mpc.ri_0,mpc.ri_k,mpc.ri_ter,...
    y_cnstr.min_limit,y_cnstr.max_limit,...
    y_cnstr.min_0,y_cnstr.min,y_cnstr.min_ter,...
    y_cnstr.max_0,y_cnstr.max,y_cnstr.max_ter,...
    y_cnstr.min_ineqRow_0,y_cnstr.min_ineqRow_k,y_cnstr.min_ineqRow_ter,...
    y_cnstr.max_ineqRow_0,y_cnstr.max_ineqRow_k,y_cnstr.max_ineqRow_ter,...
    y_cnstr.use_k0,y_cnstr.use_ter,mpc.y_0,mpc.y,mpc.y_ter);
end

function [ri_0,ri_k,ri_ter] = inequality_output_local(ri_0,ri_k,ri_ter,...
    min_limit,max_limit,min_0,min_k,min_ter,max_0,max_k,max_ter,...
    min_ineqRow_0,min_ineqRow_k,min_ineqRow_ter,...
    max_ineqRow_0,max_ineqRow_k,max_ineqRow_ter,use_k0,use_ter,y_0,y,y_ter)
if ~isempty(min_limit)
    if ~isempty(use_k0)
        ri_0(min_ineqRow_0) = ri_0(min_ineqRow_0) + min_0-y_0;
    end
    ri_k(min_ineqRow_k,:) = ri_k(min_ineqRow_k,:) + min_k-y;
    if ~isempty(use_ter)
        ri_ter(min_ineqRow_ter) = ri_ter(min_ineqRow_ter) + min_ter-y_ter;
    end
end
if ~isempty(max_limit)
    if ~isempty(use_k0)
        ri_0(max_ineqRow_0) = ri_0(max_ineqRow_0) + y_0-max_0;
    end
    ri_k(max_ineqRow_k,:) = ri_k(max_ineqRow_k,:) + y-max_k;
    if ~isempty(use_ter)
        ri_ter(max_ineqRow_ter) = ri_ter(max_ineqRow_ter) + y_ter-max_ter;
    end
end
end
