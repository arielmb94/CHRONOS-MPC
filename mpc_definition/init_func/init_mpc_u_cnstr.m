function mpc = init_mpc_u_cnstr(mpc,u_min,u_max)

mpc.u_min = u_min;
mpc.u_max = u_max;

% Control box constraints
[mpc.gradUmin,mpc.gradUmax] = genGradU(mpc.N_ctr_hor,...
                                    mpc.Nx,mpc.Nu,mpc.nx,mpc.nu);

if ~isempty(mpc.u_min)
    [mpc.hessUmin,mi] = genHessIneq(mpc.gradUmin);
    mpc.m = mpc.m+mi;
end
if ~isempty(mpc.u_max)
    [mpc.hessUmax,mi] = genHessIneq(mpc.gradUmax);
    mpc.m = mpc.m+mi;
end

end