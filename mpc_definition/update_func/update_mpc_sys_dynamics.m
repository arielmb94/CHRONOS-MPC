function mpc = update_mpc_sys_dynamics(mpc,A,B,Bd)

if ~isempty(A) && ~isempty(B)
    mpc.A = A;
    mpc.B = B;

    % A equality contraint 
    mpc.Aeq = genEqualities(A,B,mpc.N,mpc.N_ctr_hor,mpc.Nx,mpc.Nu,mpc.nx,mpc.nu);
end

if ~isempty(Bd)   
    mpc.Bd = Bd;
end

end