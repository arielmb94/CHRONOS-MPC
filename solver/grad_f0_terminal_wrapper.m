function mpc = grad_f0_terminal_wrapper(mpc)
mpc.grad_xe_f0_ter=grad_f0_terminal_local(mpc.grad_xe_f0_ter,mpc.x_col,mpc.P2,mpc.xN_ref,mpc.x_ter);
end
function grad_xe_f0_ter = grad_f0_terminal_local(grad_xe_f0_ter,x_idx,P2,xN_ref,x_ter)
grad_xe_f0_ter(x_idx) = grad_xe_f0_ter(x_idx) - P2*(xN_ref-x_ter);
end
