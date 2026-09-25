function mpc = build_optimization_variables(mpc)
% BUILD_OPTIMIZATION_VARIABLES Allocate the fixed stage-local solver workspace.

N = mpc.N;
nse = mpc.nx + double(~isempty(mpc.has_du))*mpc.nu;
ng_0 = mpc.ng_k(1);
ng_k = mpc.ng_k(2);
ng_ter = mpc.ng_k(3);
nv_0 = mpc.nv_k(1);
nv_k = mpc.nv_k(2);
nv_ter = mpc.nv_k(3);

mpc.nse = nse;
mpc.nvar = N*(nse + mpc.nu);
mpc.n = mpc.nvar + ng_0 + (N-1)*ng_k + ng_ter + ...
        nv_0 + (N-1)*nv_k + nv_ter;

mpc.s = zeros(mpc.nx,N);
mpc.s_ter = zeros(mpc.nx,1);
mpc.u = zeros(mpc.nu,N);
if ~isempty(mpc.has_du)
    mpc.su = zeros(mpc.nu,N);
    mpc.du = zeros(mpc.nu,N);
else
    mpc.su = [];
    mpc.du = [];
end

mpc.delta_u = zeros(mpc.nu,N);
mpc.delta_se = zeros(nse,N);

mpc.ru_0 = zeros(mpc.nu,1);
mpc.ru_k = zeros(mpc.nu,N-1);
mpc.rse_k = zeros(nse,N-1);
mpc.rse_ter = zeros(nse,1);

mpc.ru_hat_0 = zeros(mpc.nu,1);

mpc.ri_0 = zeros(ng_0,1);
mpc.ri_k = zeros(ng_k,N-1);
mpc.ri_ter = zeros(ng_ter,1);
mpc.iS_0 = zeros(ng_0,1);
mpc.iS_k = zeros(ng_k,N-1);
mpc.iS_ter = zeros(ng_ter,1);
mpc.iS_ri_hat_0 = zeros(ng_0,1);
mpc.iS_ri_hat_k = zeros(ng_k,N-1);
mpc.iS_ri_hat_ter = zeros(ng_ter,1);

mpc.g_0 = zeros(ng_0,1);
mpc.g_k = zeros(ng_k,N-1);
mpc.g_ter = zeros(ng_ter,1);
mpc.g2_0 = zeros(ng_0,1);
mpc.g2_k = zeros(ng_k,N-1);
mpc.g2_ter = zeros(ng_ter,1);
mpc.delta_g_0 = zeros(ng_0,1);
mpc.delta_g_k = zeros(ng_k,N-1);
mpc.delta_g_ter = zeros(ng_ter,1);

mpc.v_0 = zeros(nv_0,1);
mpc.v_k = zeros(nv_k,N-1);
mpc.v_ter = zeros(nv_ter,1);
mpc.rv_v2_0 = zeros(nv_0,1);
mpc.rv_v2_k = zeros(nv_k,N-1);
mpc.rv_v2_ter = zeros(nv_ter,1);
mpc.v2_0 = zeros(nv_0,1);
mpc.v2_k = zeros(nv_k,N-1);
mpc.v2_ter = zeros(nv_ter,1);
mpc.delta_v_0 = zeros(nv_0,1);
mpc.delta_v_k = zeros(nv_k,N-1);
mpc.delta_v_ter = zeros(nv_ter,1);

end
