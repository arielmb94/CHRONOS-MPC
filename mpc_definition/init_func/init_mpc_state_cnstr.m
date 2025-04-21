function mpc = init_mpc_state_cnstr(mpc,x_min,x_max)

mpc.x_min = x_min;
mpc.x_max = x_max;

[mpc.gradXmin,mpc.gradXmax] = genGradX(mpc.N,mpc.N_ctr_hor,...
                                mpc.Nx,mpc.Nu,mpc.nx,mpc.nu);

if ~isempty(mpc.x_min)
    [mpc.hessXmin,mi] = genHessIneq(mpc.gradXmin);
    mpc.m = mpc.m+mi;
end
if ~isempty(mpc.x_max)
    [mpc.hessXmax,mi] = genHessIneq(mpc.gradXmax);
    mpc.m = mpc.m+mi;
end

end