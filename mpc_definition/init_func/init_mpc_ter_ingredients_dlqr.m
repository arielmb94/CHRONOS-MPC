function mpc = init_mpc_ter_ingredients_dlqr(mpc,Qx,Ru,...
    terminal_constraint,x_ref_is_y)

mpc.ter_ingredients = 1;
mpc.ter_constraint = terminal_constraint;
mpc.x_ref_is_y = x_ref_is_y;

[~,P] = dlqr(mpc.A,mpc.B,Qx,Ru);

mpc.P = P;

if isempty(mpc.hessCost)
    mpc.hessCost = zeros(mpc.Nu+mpc.Nx);
end

mpc.hessTerminalCost = zeros(mpc.Nx+mpc.Nu);
mpc.hessTerminalCost(end-mpc.nx+1: end,end-mpc.nx+1 : end) = 2*P;

mpc.hessCost = mpc.hessCost + mpc.hessTerminalCost;

end