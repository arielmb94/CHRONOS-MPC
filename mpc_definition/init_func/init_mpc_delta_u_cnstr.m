function mpc = init_mpc_delta_u_cnstr(mpc,du_min,du_max)

mpc.du_min = du_min;
mpc.du_max = du_max;

% Differential Control box constraints
[mpc.gradDeltaUmin,mpc.gradDeltaUmax] = genGradDeltaU(mpc.N_ctr_hor,...
                                                mpc.Nx,mpc.Nu,mpc.nx,mpc.nu);

if ~isempty(mpc.du_min)
    [mpc.hessDeltaUmin,mi] = genHessIneq(mpc.gradDeltaUmin);
    mpc.m = mpc.m+mi;
end
if ~isempty(mpc.du_max)
    [mpc.hessDeltaUmax,mi] = genHessIneq(mpc.gradDeltaUmax);
    mpc.m = mpc.m+mi;
end

end