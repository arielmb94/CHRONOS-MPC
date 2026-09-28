function mpc = genEqualities(mpc)

nx = mpc.nx;
nu = mpc.nu;
nse = mpc.nse;

s_col = 1:nx;
if ~isempty(mpc.has_du)
    su_col = nx+1:nse;
else
    su_col = [];
end
u_col = nse+1:nse+nu;

mpc.s_col = s_col;
mpc.su_col = su_col;
mpc.se_col = [s_col su_col];
mpc.u_col = u_col;

mpc.rp_0 = zeros(nse,1);
mpc.rp_k = zeros(nse,mpc.N-1);
mpc.beq_0 = zeros(nx,1);
mpc.beq_k = zeros(nx,mpc.N-1);

end
