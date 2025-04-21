function mpc = init_mpc_ter_state_cnstr(mpc,x_ter_min,x_ter_max)

mpc.x_ter_min = x_ter_min; 
mpc.x_ter_max = x_ter_max;

[mpc.gradXtermin,mpc.gradXtermax] = genGradXter(mpc.N,mpc.N_ctr_hor,...
    mpc.Nx,mpc.Nu,mpc.nx,mpc.nu);

if ~isempty(mpc.x_ter_min)
    [mpc.hessXtermin,mi] = genHessIneq(mpc.gradXtermin);
    mpc.m = mpc.m+mi;
end
if ~isempty(mpc.x_ter_max)
    [mpc.hessXtermax,mi] = genHessIneq(mpc.gradXtermax);
    mpc.m = mpc.m+mi;
end

end