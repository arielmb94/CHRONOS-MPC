function mpc = init_mpc_system(mpc,A,B,Bd,C,D,Dd)

mpc.A = A;
mpc.B = B;
mpc.Bd = Bd;
mpc.C = C;
mpc.D = D;
mpc.Dd = Dd;

mpc.nx = size(A,1);  %number of states
mpc.nu = size(B,2);  %number of control inputs
mpc.nd = size(Bd,2);  %number of disturbance inputs
mpc.ny = size(C,1);  %number of measurements

mpc.Nx = mpc.N*mpc.nx;
mpc.Nu = mpc.N_ctr_hor*mpc.nu;
mpc.Nd = mpc.N*mpc.nd;

if mpc.D == 0
    mpc.Ny = (mpc.N-1)*mpc.ny;
else
    mpc.Ny = mpc.N*mpc.ny;
end

% A equality contraint (b equality constraints depends on x0 and d(k)
mpc.Aeq = genEqualities(A,B,mpc.N,mpc.N_ctr_hor,mpc.Nx,mpc.Nu,mpc.nx,mpc.nu);
mpc.beq = zeros(size(mpc.Aeq,1),1);

end