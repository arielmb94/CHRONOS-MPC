function mpc = genEqualities(mpc)

nx = mpc.nx;
nu = mpc.nu;
nxe = mpc.nxe;

x_col = 1:nx;
if ~isempty(mpc.has_du)
    xu_col = nx+1:nxe;
else
    xu_col = [];
end
u_col = nxe+1:nxe+nu;

mpc.x_col = x_col;
mpc.xu_col = xu_col;
mpc.xe_col = [x_col xu_col];
mpc.u_col = u_col;

mpc.rp_0 = zeros(nxe,1);
mpc.rp_k = zeros(nxe,mpc.N-1);
mpc.beq_0 = zeros(nx,1);
mpc.beq_k = zeros(nx,mpc.N-1);

end
